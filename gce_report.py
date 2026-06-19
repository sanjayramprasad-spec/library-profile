"""GCE clone report - protein-level analysis of evolved/library clones.

For genetic code expansion (GCE) work: takes a batch of whole-plasmid results
(e.g. Plasmidsaurus) for clones from a synthetase/tRNA library or evolution
campaign, and reports what matters at the PROTEIN level - not just nucleotides.

It builds on compare_plasmid.py (imported, not re-implemented): each clone is
aligned to the parent map, then this script adds four GCE-specific layers:

  1. PROTEIN VARIANTS   - amino-acid changes in the gene of interest vs the
                          parent (e.g. Y32G, L65V), plus changes OUTSIDE the
                          designed positions (often where activity comes from).
  2. AMBER / STOP AUDIT - finds every stop codon in the gene, marks the designed
                          amber (TAG) sites as expected, and flags unintended
                          stops - distinguishing suppressible TAG from hard
                          TAA/TGA. Also flags a configured amber that is absent
                          (e.g. a numbering offset, or a reverted site).
  3. NNK / LIBRARY QC   - at each randomized position: the observed residue,
                          codon, and whether the codon obeys the scheme (NNK).
  4. CAMPAIGN ROLL-UP   - across reliable clones: a clone x mutation matrix, unique
                          genotypes, and CONVERGENT mutations. Recurrence is a
                          HYPOTHESIS (a lead), not proof: it implies selection only
                          across INDEPENDENT isolates, and sequence convergence is
                          not fitness convergence - confirm with activity data.

USAGE
    python gce_report.py PARENT.dna  RESULTS_FOLDER  --config campaign.txt
    python gce_report.py PARENT.dna  clone1.fasta clone2.fasta --config campaign.txt

CONFIG FILE (plain text; '#' starts a comment; one 'key: value' per line)
    gene: TargetGene           # label of the CDS in the map to analyse
    randomized: 15, 25, 35     # 1-based residue numbers randomized (NNK). optional
    amber: 40                  # designed amber (TAG) sites, residue numbers. optional
    scheme: NNK                # randomization scheme. optional, default NNK

ASSUMPTION - the stop audit assumes AMBER (TAG) suppression, the standard for the
    vast majority of GCE. Only TAG is treated as suppressible; TGA (opal) and TAA
    (ochre) are reported as hard stops. If your system suppresses opal/ochre,
    read those "HARD STOP" findings with that in mind.

OUTPUTS (in "gce_report/" next to the parent .dna)
    gce_variants.csv      one row per clone: mutations, pocket genotype, verdict
    mutation_matrix.csv   clones x mutated positions
    gce_summary.txt       human-readable report incl. convergent mutations
    mutation_matrix.png   heatmap of clones x positions

Dependencies: biopython, matplotlib, and compare_plasmid.py in the same folder.
"""
from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

try:
    from Bio.Seq import Seq
    # Reuse the validated alignment engine - do NOT parse its text output.
    from compare_plasmid import (read_reference, read_query, compare, _gather_queries,
                                 _flabel, warn_if_untested_biopython)
except ModuleNotFoundError as _exc:   # friendly first-run message for non-programmers
    raise SystemExit(
        f"ERROR: missing required package '{_exc.name}'.\n"
        f"  Install dependencies once with:  pip install -r requirements.txt\n"
        f"  (or:  pip install biopython matplotlib)"
    )

STOP_CODONS = {"TAA", "TAG", "TGA"}
# Which 3rd-base nucleotides a scheme allows (used to check library codons).
_SCHEME_WOBBLE = {"NNK": set("GT"), "NNS": set("GC"), "NNN": set("ACGT")}


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Campaign:
    gene: str
    randomized: tuple[int, ...]
    amber: tuple[int, ...]
    scheme: str


def _int_list(value: str) -> tuple[int, ...]:
    out = []
    for tok in value.replace(",", " ").split():
        try:
            out.append(int(tok))
        except ValueError:
            pass
    return tuple(out)


def parse_config(path: Path) -> Campaign:
    cfg = {"gene": None, "randomized": (), "amber": (), "scheme": "NNK"}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or ":" not in line:
            continue
        key, value = line.split(":", 1)
        key, value = key.strip().lower(), value.strip()
        if key == "gene":
            cfg["gene"] = value
        elif key == "randomized":
            cfg["randomized"] = _int_list(value)
        elif key == "amber":
            cfg["amber"] = _int_list(value)
        elif key == "scheme":
            cfg["scheme"] = value.upper()
    if not cfg["gene"]:
        raise ValueError("config must specify 'gene:' (the CDS label to analyse)")
    return Campaign(cfg["gene"], cfg["randomized"], cfg["amber"], cfg["scheme"])


def find_cds(ref, gene_label: str):
    """Find the CDS feature whose label matches the configured gene."""
    cds = [f for f in ref.features if f.type == "CDS"]
    for f in cds:                                   # exact label match first
        if _flabel(f).lower() == gene_label.lower():
            return f
    for f in cds:                                   # then substring
        if gene_label.lower() in _flabel(f).lower():
            return f
    return None


# --------------------------------------------------------------------------- #
# Per-clone analysis
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GeneVariant:
    resnum: int
    parent_aa: str
    variant_aa: str

    def __str__(self) -> str:                       # e.g. "Y32G", or "*187Q"
        return f"{self.parent_aa}{self.resnum}{self.variant_aa}"


@dataclass(frozen=True)
class PocketPosition:
    resnum: int
    parent_aa: str
    variant_aa: str
    codon: str
    scheme_ok: bool


@dataclass(frozen=True)
class CloneResult:
    name: str
    identity: float
    mutations: tuple[GeneVariant, ...]
    pocket: tuple[PocketPosition, ...]
    amber_findings: tuple[str, ...]
    has_cds_indel: bool
    verdict: str


def _coding(cds, genome: Seq) -> Seq:
    """Strand-aware coding sequence for a CDS, uppercased."""
    return Seq(str(cds.extract(genome)).upper())


def _apply_substitutions(ref, result) -> tuple[Seq, list]:
    """Build the variant genome by applying point mutations; collect indels.

    GCE library hits are overwhelmingly substitutions (NNK). Indels are
    returned separately so a frameshift can be flagged rather than mistranslated.
    """
    seq = list(str(ref.seq).upper())
    indels = []
    for d in result.differences:
        if d.kind == "mismatch":
            seq[d.ref_pos - 1] = d.query_base.upper()
        else:
            indels.append(d)
    return Seq("".join(seq)), indels


def _amber_audit(var_coding: Seq, var_prot: str, amber_sites: tuple[int, ...]) -> tuple[str, ...]:
    """Classify every stop in the variant ORF against the designed amber sites."""
    stops = [i + 1 for i, a in enumerate(var_prot) if a == "*"]
    configured = set(amber_sites)
    findings = []

    for site in amber_sites:                        # is each designed amber present?
        if site in stops:
            codon = str(var_coding[(site - 1) * 3:site * 3]).upper()
            tag = "TAG (amber)" if codon == "TAG" else f"{codon}"
            findings.append(f"designed amber at residue {site} present: {tag}")
        else:
            aa = var_prot[site - 1] if site - 1 < len(var_prot) else "?"
            near = [p for p in stops if abs(p - site) <= 2]
            hint = (f"  [but a stop IS present at residue {near[0]} - check numbering offset]"
                    if near else "")
            findings.append(f"designed amber at residue {site} ABSENT (residue is {aa}){hint}")

    for p in stops:                                 # any unconfigured stop?
        if p not in configured:
            codon = str(var_coding[(p - 1) * 3:p * 3]).upper()
            kind = "TAG amber (suppressible)" if codon == "TAG" else f"{codon} HARD STOP (not suppressible)"
            findings.append(f"unexpected stop at residue {p}: {kind}")
    return tuple(findings)


def _make_verdict(muts, amber_findings, has_cds_indel) -> str:
    if has_cds_indel:
        return "FRAMESHIFT - CDS indel, protein call unreliable"
    if any("HARD STOP" in f for f in amber_findings):
        return "PREMATURE HARD STOP - truncated even in a suppressor host"
    amber_issue = any("ABSENT" in f for f in amber_findings)
    parts = []
    if muts:
        parts.append(f"{len(muts)} aa change(s)")
    if amber_issue:
        parts.append("amber-site issue")
    if not parts:
        return "parent (no coding changes; designed amber intact)"
    return ", ".join(parts)


def analyze_clone(ref, query, cds, camp: Campaign) -> CloneResult:
    result = compare(ref, query)
    variant_genome, indels = _apply_substitutions(ref, result)

    parent_prot = str(_coding(cds, ref.seq).translate())
    var_coding = _coding(cds, variant_genome)
    var_prot = str(var_coding.translate())

    # 1. protein variants in the gene of interest
    muts = tuple(GeneVariant(i + 1, p, v)
                 for i, (p, v) in enumerate(zip(parent_prot, var_prot)) if p != v)

    # 3. NNK / library pocket genotype at the randomized positions
    wobble = _SCHEME_WOBBLE.get(camp.scheme, set("ACGT"))
    pocket = []
    for pos in camp.randomized:
        if pos * 3 <= len(var_coding):
            codon = str(var_coding[(pos - 1) * 3:pos * 3]).upper()
            aa = str(Seq(codon).translate()) if len(codon) == 3 else "?"
            paa = parent_prot[pos - 1] if pos - 1 < len(parent_prot) else "?"
            pocket.append(PocketPosition(pos, paa, aa, codon, codon[2] in wobble))

    # 2. amber / stop audit
    amber_findings = _amber_audit(var_coding, var_prot, camp.amber)

    # frameshift flag: any indel inside the CDS span
    cds_start, cds_end = int(cds.location.start), int(cds.location.end)
    has_cds_indel = any(cds_start < d.ref_pos <= cds_end + 1 for d in indels)

    return CloneResult(
        name=result.query_name,
        identity=result.percent_identity,
        mutations=muts,
        pocket=tuple(pocket),
        amber_findings=amber_findings,
        has_cds_indel=has_cds_indel,
        verdict=_make_verdict(muts, amber_findings, has_cds_indel),
    )


# --------------------------------------------------------------------------- #
# Campaign roll-up (across clones)
# --------------------------------------------------------------------------- #
def _reliable(clones: list[CloneResult]) -> list[CloneResult]:
    """Clones whose protein call is trustworthy. Frameshifted clones are excluded:
    their in-frame substitution list is an artifact (ISSUE-5 flagged it UNRELIABLE),
    so it must not feed the convergence/genotype roll-up (CONCERN-2)."""
    return [c for c in clones if not c.has_cds_indel]


def convergent_mutations(clones: list[CloneResult]) -> list[tuple[str, int]]:
    """Mutations ranked by how many RELIABLE clones carry them.

    Recurrence is a HYPOTHESIS, not proof of function: it implies selection only
    across INDEPENDENT isolates (siblings from one outgrowth inflate it), and
    sequence convergence is not fitness convergence - confirm with activity data.
    """
    freq = Counter(str(m) for c in _reliable(clones) for m in c.mutations)
    return sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))


def unique_genotypes(clones: list[CloneResult]) -> dict[tuple, list[str]]:
    geno: dict[tuple, list[str]] = {}
    for c in _reliable(clones):
        key = tuple(sorted(str(m) for m in c.mutations))
        geno.setdefault(key, []).append(c.name)
    return geno


def _pocket_str(c: CloneResult) -> str:
    return " ".join(f"{p.resnum}:{p.variant_aa}{'' if p.scheme_ok else '(!)'}" for p in c.pocket)


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def print_clone(c: CloneResult) -> None:
    muts = ", ".join(str(m) for m in c.mutations) or "(none)"
    print(f"\n  {c.name}   [{c.identity:.2f}% identity]")
    print(f"    verdict   : {c.verdict}")
    print(f"    mutations : {muts}")
    if c.pocket:
        print(f"    pocket    : {_pocket_str(c)}    [(!) = codon off-scheme]")
    for f in c.amber_findings:
        print(f"    amber     : {f}")


def write_outputs(clones: list[CloneResult], camp: Campaign, ref_name: str, out_dir: Path) -> None:
    import csv

    # gce_variants.csv - one row per clone
    with (out_dir / "gce_variants.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["clone", "identity_%", "verdict", "gene_mutations", "pocket_genotype",
                    "amber_status"])
        for c in clones:
            # A CDS frameshift makes the in-frame protein call meaningless past the
            # indel; flag the per-residue fields so the unreliability travels with
            # the data (a column-only consumer would otherwise miss the verdict).
            muts = "UNRELIABLE (frameshift)" if c.has_cds_indel else ";".join(str(m) for m in c.mutations)
            pocket = "UNRELIABLE (frameshift)" if c.has_cds_indel else _pocket_str(c)
            w.writerow([c.name, f"{c.identity:.2f}", c.verdict, muts, pocket,
                        " | ".join(c.amber_findings)])

    # mutation_matrix.csv - clones x mutated positions (columns from RELIABLE clones)
    positions = sorted({m.resnum for c in clones if not c.has_cds_indel for m in c.mutations}
                       | set(camp.randomized))
    with (out_dir / "mutation_matrix.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["clone"] + [str(p) for p in positions])
        for c in clones:
            if c.has_cds_indel:                     # frameshift: residue calls unreliable
                w.writerow([c.name] + ["fs"] * len(positions))
            else:
                at = {m.resnum: m.variant_aa for m in c.mutations}
                w.writerow([c.name] + [at.get(p, ".") for p in positions])

    # gce_summary.txt
    lines = [f"GCE campaign report  |  parent: {ref_name}  |  gene: {camp.gene}",
             f"clones analysed: {len(clones)}", "=" * 70, ""]
    rel = _reliable(clones)
    n_excl = len(clones) - len(rel)
    excl_note = f"; {n_excl} frameshift clone(s) excluded from roll-up" if n_excl else ""
    geno = unique_genotypes(clones)
    lines.append(f"UNIQUE GENOTYPES: {len(geno)} distinct (from {len(rel)} clones{excl_note})")
    for key, names in sorted(geno.items(), key=lambda kv: -len(kv[1])):
        muts = ", ".join(key) if key else "(parent / no changes)"
        lines.append(f"  x{len(names):<3} {muts}    [{', '.join(names)}]")
    lines += ["", "CONVERGENT MUTATIONS - a HYPOTHESIS, not proof (recurrence implies selection",
              "only across INDEPENDENT isolates; sequence convergence is not fitness - confirm",
              "with activity data):"]
    conv = convergent_mutations(clones)
    if conv:
        for mut, n in conv:
            star = "  <- convergent" if n > 1 else ""
            lines.append(f"  {mut:<10} in {n}/{len(rel)} clones{star}")
    else:
        lines.append("  (none)")
    lines += ["", "PER-CLONE:"]
    for c in clones:
        lines.append(f"  {c.name}: {c.verdict}")
        if c.amber_findings:
            for f in c.amber_findings:
                lines.append(f"       amber: {f}")
    (out_dir / "gce_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_matrix(clones: list[CloneResult], camp: Campaign, out_path: Path) -> None:
    """Heatmap: clones (rows) x positions (cols); cells show the variant residue."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    positions = sorted({m.resnum for c in clones if not c.has_cds_indel for m in c.mutations}
                       | set(camp.randomized))
    if not positions or not clones:
        return
    fig, ax = plt.subplots(figsize=(1.2 + 0.55 * len(positions), 1.0 + 0.45 * len(clones)))
    for r, c in enumerate(clones):
        if c.has_cds_indel:                         # frameshift: in-frame residues unreliable
            for col in range(len(positions)):
                ax.add_patch(plt.Rectangle((col, r), 1, 1, facecolor="#999999", edgecolor="white"))
            ax.text(len(positions) / 2, r + 0.5, "frameshift - call unreliable",
                    ha="center", va="center", fontsize=7, color="white", style="italic")
            continue
        at = {m.resnum: m.variant_aa for m in c.mutations}
        for col, pos in enumerate(positions):
            aa = at.get(pos)
            changed = aa is not None
            ax.add_patch(plt.Rectangle((col, r), 1, 1, facecolor=("#d62728" if changed else "#f0f0f0"),
                                       edgecolor="white"))
            if changed:
                ax.text(col + 0.5, r + 0.5, aa, ha="center", va="center", fontsize=8,
                        color="white", weight="bold")
    ax.set_xlim(0, len(positions))
    ax.set_ylim(0, len(clones))
    ax.set_xticks([i + 0.5 for i in range(len(positions))])
    ax.set_xticklabels([str(p) for p in positions], fontsize=8)
    ax.set_yticks([i + 0.5 for i in range(len(clones))])
    ax.set_yticklabels([c.name[:28] for c in clones], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"residue position in {camp.gene}")
    ax.set_title(f"{camp.gene}: clone x mutation map ({len(clones)} clones)")
    ax.tick_params(length=0)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _parse_args(argv: list[str]) -> tuple[list[str], Path | None]:
    positional, config = [], None
    i = 0
    while i < len(argv):
        if argv[i] == "--config":
            config = Path(argv[i + 1]) if i + 1 < len(argv) else None
            i += 2
        else:
            positional.append(argv[i])
            i += 1
    return positional, config


def main(argv: list[str]) -> int:
    positional, config_path = _parse_args(argv)
    if len(positional) < 2 or config_path is None:
        print(__doc__)
        return 1
    warn_if_untested_biopython()
    if not config_path.exists():
        print(f"ERROR: config file not found: {config_path}")
        return 1

    ref_path = Path(positional[0])
    if not ref_path.exists():
        print(f"ERROR: parent map not found: {ref_path}")
        return 1

    camp = parse_config(config_path)
    ref = read_reference(ref_path)
    cds = find_cds(ref, camp.gene)
    if cds is None:
        labels = sorted({_flabel(f) for f in ref.features if f.type == "CDS"})
        print(f"ERROR: no CDS labelled '{camp.gene}' in the map. CDS labels present: {labels}")
        return 1
    # Translation assumes the CDS starts on a full codon (codon_start=1). A partial
    # CDS (codon_start != 1) would shift every residue number and call; warn loudly
    # rather than report wrong protein positions. (SnapGene maps always use 1.)
    if int(cds.qualifiers.get("codon_start", [1])[0]) != 1:
        print(f"  WARNING: CDS '{camp.gene}' has codon_start != 1 (partial/shifted frame); "
              f"residue numbers and protein calls assume codon_start=1 and may be WRONG.")

    queries = _gather_queries(positional[1:])
    if not queries:
        print("ERROR: no clone sequences found. Pass .fasta/.gbk/.ab1 files or a folder.")
        return 1

    out_dir = ref_path.parent / "gce_report"
    out_dir.mkdir(exist_ok=True)
    print(f"Parent: {ref.name}  |  gene: {camp.gene} ({len(_coding(cds, ref.seq)) // 3} aa)")
    print(f"Randomized: {list(camp.randomized) or '-'}   Amber: {list(camp.amber) or '-'}   "
          f"Scheme: {camp.scheme}")
    print(f"Analysing {len(queries)} clone(s); outputs -> {out_dir}")

    clones = []
    for qpath in queries:
        try:
            c = analyze_clone(ref, read_query(qpath), cds, camp)
        except Exception as e:               # one bad clone must not stop the batch
            print(f"\n  SKIPPED {qpath.name}: {e}")
            continue
        print_clone(c)
        clones.append(c)

    if not clones:
        print("\nNo clones analysed.")
        return 1

    write_outputs(clones, camp, ref.name, out_dir)
    plot_matrix(clones, camp, out_dir / "mutation_matrix.png")

    rel = _reliable(clones)
    n_excl = len(clones) - len(rel)
    print(f"\n{'='*70}")
    print("CONVERGENT MUTATIONS - a HYPOTHESIS (selection only across INDEPENDENT isolates;")
    print("  sequence convergence is not fitness - confirm with activity data):")
    for mut, n in convergent_mutations(clones):
        print(f"  {mut:<10} in {n}/{len(rel)} clones" + ("   <- convergent" if n > 1 else ""))
    excl = f"  ({n_excl} frameshift clone(s) excluded)" if n_excl else ""
    print(f"\nUnique genotypes: {len(unique_genotypes(clones))} from {len(rel)} clones{excl}")
    print(f"Outputs -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
