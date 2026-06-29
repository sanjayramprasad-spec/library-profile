"""Profile a CHIMERA LIBRARY from pooled long reads — which source contributed
each domain of every construct, and how the combinations are distributed.

Third companion to compare_plasmid.py / cloning_report.py. Those validate ONE
construct against ONE map. This one is for a *combinatorial library*: many
constructs, each assembled from interchangeable domain fragments taken from a
panel of source genes (e.g. an N-terminal, a middle, and a C-terminal domain
each drawn from a different source and fused at conserved junctions). One pooled
Nanopore run (Plasmidsaurus) gives thousands of reads; per molecule we want the
*domain mosaic* — e.g. N=srcA | middle=srcC | C=srcB — then the library-wide
composition, abundance skew, and dropouts.

WHY NO ALIGNER (the design choice)
    A read of a chimera matches NO single reference end-to-end, so a "best hit"
    is meaningless. But the source alleles within each domain are 20-45% divergent
    from one another, so each allele carries hundreds of k-mers unique to it.
    Counting which allele's private k-mers a read contains, per domain, classifies
    every domain unambiguously even at ~5% Nanopore error — no alignment, and it
    stays pure-Python (biopython only), so it installs and runs on Windows where
    minimap2/mappy do not build. The references' shared backbone and conserved
    junctions need no special handling: backbone k-mers are common to all alleles
    (never private), so they simply never vote.

HOW IT WORKS (the pipeline, in order)
    1. read the multi-FASTA of full-length source references (one per source) and
       align them to each other (center-star) -> reference_msa().
    2. DETECT DOMAINS automatically: long perfectly-conserved runs are the
       assembly junctions; the variable blocks between them are the domains
       (their count is discovered, not assumed) -> detect_domains().
    3. cut each reference into its per-domain alleles -> extract_alleles().
    4. for each domain, collect the k-mers PRIVATE to one allele (the markers)
       -> build_markers(). Backbone / conserved-junction k-mers are shared, so
       they are dropped automatically and cannot vote.
    5. for each read, in whichever orientation fits, tally private-marker hits
       per allele per domain; the top allele wins that domain -> call_read().
    6. HONESTY GATE: a domain is only called when its winner clears an absolute
       marker floor AND beats the runner-up by a margin; otherwise it is
       'unassigned' (a partial read, or a source NOT in your reference panel).
       A genotype counts toward composition only if every domain is called.
    7. report: per-read TSV, composition + per-domain usage TSV, a domain x source
       usage heatmap and a top-genotypes bar, and a QC summary -> main().

USAGE
    python library_profile.py REFERENCES.fasta READS.fastq [READS2.fastq ...]
    python library_profile.py REFERENCES.fasta a_folder_of_fastqs
    options:
      --k 15                 k-mer length (odd; larger = stricter, fewer hits)
      --min-markers 10       absolute private-marker floor to call a domain
      --margin 3.0           winner must beat runner-up by this ratio (by marker fraction)
      --names N,cat,C        friendly domain names (else dom1, dom2, ... in order)
      --anchor-min 20        min length of a conserved run treated as a junction
      --expected FILE.tsv    optional designed-combination list for a coverage report

OUTPUTS (in "library_profile/" next to the references file)
    per_read.tsv            one row per read: orientation, per-domain call, status
    composition.tsv         each observed source-combination and its read share
    domain_usage.tsv        per domain, how often each source allele was used
    domain_usage.png        domain x source usage heatmap
    top_genotypes.png       the most abundant combinations
    library_summary.txt     honesty-gated QC: yields, dropouts, skew, caveats

Dependencies: biopython, matplotlib   (see requirements.txt)
"""
from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
except ModuleNotFoundError as _exc:   # friendly first-run message for non-programmers
    raise SystemExit(
        f"ERROR: missing required package '{_exc.name}'.\n"
        f"  Install dependencies once with:  pip install -r requirements.txt\n"
        f"  (or:  pip install biopython matplotlib)"
    )

# --- Defaults (all overridable on the command line) ------------------------ #
DEFAULT_K = 15            # k-mer length. 15 is long enough to be source-specific,
                          # short enough that ~5% Nanopore error still leaves many
                          # error-free windows per domain (measured: 300-600/allele).
DEFAULT_MIN_MARKERS = 10  # below this many private-marker hits a domain is 'unassigned'
DEFAULT_MARGIN = 3.0      # winner's marker-fraction must be >= this x the runner-up's
DEFAULT_ANCHOR_MIN = 20   # a conserved run >= this long is treated as an assembly junction
MIN_DOMAIN_BP = 30        # ignore variable specks shorter than this when splitting domains

_READ_FORMATS = {".fastq": "fastq", ".fq": "fastq", ".fasta": "fasta",
                 ".fa": "fasta", ".fna": "fasta", ".gz": None}
_COMPL = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _rc(seq: str) -> str:
    return seq.translate(_COMPL)[::-1]


# Read status, worst-wins. 'partial' outranks 'ambiguous' because a missing /
# unknown-source domain is the more actionable signal (add the missing source).
_STATUS_RANK = {"complete": 0, "ambiguous": 1, "partial": 2}


def _worse_status(current: str, candidate: str) -> str:
    return candidate if _STATUS_RANK[candidate] > _STATUS_RANK[current] else current


# --------------------------------------------------------------------------- #
# Reference parsing, alignment, and domain detection
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Domain:
    """A variable block (one domain) in center-reference coordinates."""
    name: str
    start: int       # inclusive, center coords
    end: int         # exclusive, center coords


@dataclass(frozen=True)
class ReferencePanel:
    """The parsed library design: source names, per-domain alleles, and the
    private-marker index used to classify reads."""
    sources: tuple[str, ...]
    domains: tuple[Domain, ...]
    alleles: dict                      # domain_name -> {source -> allele_seq}
    marker_pool: dict                  # domain_name -> {source -> n_private_kmers}
    index: dict                        # kmer -> (domain_name, source)
    k: int


def read_references(path: Path) -> list:
    """Read the multi-FASTA of full-length source references (>=2 records)."""
    recs = list(SeqIO.parse(str(path), "fasta"))
    if len(recs) < 2:
        raise ValueError(
            f"need >=2 source references in '{path.name}', found {len(recs)}. "
            "This tool profiles a library built from a PANEL of sources."
        )
    seen: set[str] = set()
    for r in recs:
        if r.id in seen:
            raise ValueError(f"duplicate reference name '{r.id}' in {path.name}")
        seen.add(r.id)
    return recs


def _make_aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -5.0
    a.extend_gap_score = -0.5
    return a


def reference_msa(refs: list) -> tuple[str, dict]:
    """Center-star alignment: align every reference to the longest one and read
    each off in the center's coordinates.

    Returns (center_id, cols) where cols[i] maps each source -> the base aligned
    to center position i (or '-' where that source has a deletion). Insertions
    relative to the center are dropped: domain identity comes from substitutions
    within a domain's core, so center coordinates are a sufficient common frame.
    """
    seqs = {r.id: str(r.seq).upper() for r in refs}
    center = max(seqs, key=lambda i: len(seqs[i]))
    c = seqs[center]
    cols: list[dict] = [{center: ch} for ch in c]
    aligner = _make_aligner()
    for i, s in seqs.items():
        if i == center:
            continue
        aln = aligner.align(c, s)[0]
        # str() the aligned rows: across Biopython versions indexing an Alignment
        # row has returned either a string or an array; str() is stable (the
        # sibling compare_plasmid.py relies on the same cast).
        a_center, a_other = str(aln[0]), str(aln[1])
        ci = -1
        for ac, bc in zip(a_center, a_other):
            if ac != "-":
                ci += 1
                cols[ci][i] = bc        # '-' where the other source is deleted here
    return center, {"cols": cols, "len": len(c), "center_seq": c}


def _conserved_mask(msa: dict, sources: tuple[str, ...]) -> list[bool]:
    """True at center columns where every source carries the same base (no gap)."""
    cols = msa["cols"]
    mask = []
    for col in cols:
        bases = {col.get(s) for s in sources}
        mask.append(len(bases) == 1 and None not in bases and "-" not in bases)
    return mask


def detect_domains(msa: dict, sources: tuple[str, ...], *,
                   anchor_min: int = DEFAULT_ANCHOR_MIN,
                   min_domain_bp: int = MIN_DOMAIN_BP,
                   names: tuple[str, ...] | None = None) -> tuple[Domain, ...]:
    """Discover domains as the variable blocks between conserved assembly junctions.

    The variable region is the span from the first to the last non-conserved
    column (everything outside it is shared backbone). Inside it, any
    perfectly-conserved run >= anchor_min is a junction; the variable stretches
    it separates are the domains. The number of domains is discovered, not fixed.
    """
    mask = _conserved_mask(msa, sources)
    var_positions = [i for i, cons in enumerate(mask) if not cons]
    if not var_positions:
        raise ValueError("references are identical across their whole length; "
                         "nothing to profile (no variable domains).")
    v0, v1 = var_positions[0], var_positions[-1] + 1

    # Conserved runs (length >= anchor_min) strictly inside the variable window
    # are the junctions that cut it into domains.
    cut_points: list[tuple[int, int]] = []
    run_start = None
    for i in range(v0, v1 + 1):
        cons = i < v1 and mask[i]
        if cons and run_start is None:
            run_start = i
        elif not cons and run_start is not None:
            if i - run_start >= anchor_min:
                cut_points.append((run_start, i))
            run_start = None

    # Build domain spans = variable window minus the junction cuts.
    spans: list[tuple[int, int]] = []
    cursor = v0
    for js, je in cut_points:
        if js - cursor >= min_domain_bp:
            spans.append((cursor, js))
        cursor = je
    if v1 - cursor >= min_domain_bp:
        spans.append((cursor, v1))

    if names is not None and len(names) != len(spans):
        raise ValueError(
            f"--names gave {len(names)} names but {len(spans)} domains were "
            f"detected: {[(s, e) for s, e in spans]}. Omit --names or match the count."
        )
    labels = names if names is not None else tuple(f"dom{i+1}" for i in range(len(spans)))
    return tuple(Domain(lbl, s, e) for lbl, (s, e) in zip(labels, spans))


def extract_alleles(msa: dict, sources: tuple[str, ...],
                    domains: tuple[Domain, ...]) -> dict:
    """Cut each source reference into its per-domain allele (gaps removed)."""
    cols = msa["cols"]
    alleles: dict = {d.name: {} for d in domains}
    for d in domains:
        for s in sources:
            bases = [cols[i].get(s, "-") for i in range(d.start, d.end)]
            alleles[d.name][s] = "".join(b for b in bases if b not in ("-", None))
    return alleles


# --------------------------------------------------------------------------- #
# Marker construction and read classification
# --------------------------------------------------------------------------- #
def _kmers(seq: str, k: int) -> set:
    return {km for i in range(len(seq) - k + 1)
            if "N" not in (km := seq[i:i + k])}


def build_markers(alleles: dict, domains: tuple[Domain, ...], sources: tuple[str, ...],
                  k: int) -> tuple[dict, dict]:
    """Index k-mers that are PRIVATE to one (domain, source).

    A marker must be unique to a single source within its domain AND not occur as
    a marker for any other (domain, source) anywhere in the panel. Shared backbone
    and conserved-junction k-mers therefore drop out and never vote.

    Returns (index, pool) where index[kmer] = (domain_name, source) and
    pool[domain_name][source] = how many private markers that allele has.
    """
    # First pass: per domain, which sources carry each k-mer.
    owners: dict = {}            # kmer -> set of (domain_name, source)
    for d in domains:
        seen_in_domain: dict = {}
        for s in sources:
            for km in _kmers(alleles[d.name][s], k):
                seen_in_domain.setdefault(km, set()).add(s)
        for km, srcs in seen_in_domain.items():
            if len(srcs) == 1:                       # private within this domain
                (only,) = tuple(srcs)
                owners.setdefault(km, set()).add((d.name, only))

    # Keep only k-mers owned by exactly one (domain, source) across the whole panel.
    index: dict = {}
    pool: dict = {d.name: {s: 0 for s in sources} for d in domains}
    for km, who in owners.items():
        if len(who) == 1:
            (dom, src), = tuple(who)
            index[km] = (dom, src)
            pool[dom][src] += 1
    return index, pool


def build_panel(refs: list, *, k: int = DEFAULT_K,
                anchor_min: int = DEFAULT_ANCHOR_MIN,
                names: tuple[str, ...] | None = None) -> ReferencePanel:
    """Full reference-side setup: align, detect domains, cut alleles, build markers."""
    if k % 2 == 0:
        raise ValueError("k must be odd so a k-mer and its reverse complement differ.")
    sources = tuple(r.id for r in refs)
    _, msa = reference_msa(refs)
    domains = detect_domains(msa, sources, anchor_min=anchor_min, names=names)
    alleles = extract_alleles(msa, sources, domains)
    index, pool = build_markers(alleles, domains, sources, k)
    thin = [d for d in domains if max(pool[d.name].values()) < 1]
    if thin:
        raise ValueError(
            f"domain(s) {[d.name for d in thin]} have no source-specific {k}-mers; "
            "the sources may be too similar there to tell apart."
        )
    return ReferencePanel(sources, domains, alleles, pool, index, k)


@dataclass(frozen=True)
class DomainCall:
    source: str | None       # winning source, or None if unassigned
    hits: int                # private-marker hits for the winner
    fraction: float          # hits / winner's marker-pool size (0..1, a soft identity)
    runner_up: str | None
    runner_fraction: float


@dataclass(frozen=True)
class ReadCall:
    read_id: str
    length: int
    strand: str              # '+', '-', or '?' (no markers either way)
    calls: dict              # domain_name -> DomainCall
    status: str              # 'complete' | 'partial' | 'ambiguous'

    @property
    def genotype(self) -> tuple | None:
        if self.status != "complete":
            return None
        return tuple(self.calls[d].source for d in self.calls)


def _tally(read_kmers: set, panel: ReferencePanel) -> dict:
    """One pass over a read's k-mers -> per-domain per-source hit counts."""
    counts = {d.name: Counter() for d in panel.domains}
    for km in read_kmers:
        hit = panel.index.get(km)
        if hit is not None:
            dom, src = hit
            counts[dom][src] += 1
    return counts


def call_read(read_id: str, seq: str, panel: ReferencePanel, *,
              min_markers: int = DEFAULT_MIN_MARKERS,
              margin: float = DEFAULT_MARGIN) -> ReadCall:
    """Classify one read into a per-domain source mosaic, with honesty gating."""
    seq = seq.upper()
    fwd = _kmers(seq, panel.k)
    rev = _kmers(_rc(seq), panel.k)
    fwd_counts, rev_counts = _tally(fwd, panel), _tally(rev, panel)
    fwd_total = sum(sum(c.values()) for c in fwd_counts.values())
    rev_total = sum(sum(c.values()) for c in rev_counts.values())
    if fwd_total == 0 and rev_total == 0:
        strand, counts = "?", fwd_counts
    elif fwd_total >= rev_total:
        strand, counts = "+", fwd_counts
    else:
        strand, counts = "-", rev_counts

    calls: dict = {}
    status = "complete"
    for d in panel.domains:
        c = counts[d.name]
        pool = panel.marker_pool[d.name]
        ranked = sorted(c.items(), key=lambda kv: kv[1], reverse=True)
        best_src, best_hits = (ranked[0] if ranked else (None, 0))
        best_frac = best_hits / pool[best_src] if best_src and pool[best_src] else 0.0
        run_src, run_hits = (ranked[1] if len(ranked) > 1 else (None, 0))
        run_frac = run_hits / pool[run_src] if run_src and pool[run_src] else 0.0
        # A runner-up only counts toward the margin if it has enough private
        # markers to be trustworthy; a source with a tiny marker pool can score
        # run_frac=1.0 from a few chance hits and wrongly veto a clear winner.
        run_solid = bool(run_src) and pool[run_src] >= min_markers

        if best_hits < min_markers:
            calls[d.name] = DomainCall(None, best_hits, best_frac, run_src, run_frac)
            status = _worse_status(status, "partial")
        elif run_solid and run_frac > 0 and best_frac < margin * run_frac:
            calls[d.name] = DomainCall(None, best_hits, best_frac, run_src, run_frac)
            status = _worse_status(status, "ambiguous")
        else:
            calls[d.name] = DomainCall(best_src, best_hits, best_frac, run_src, run_frac)
    return ReadCall(read_id, len(seq), strand, calls, status)


# --------------------------------------------------------------------------- #
# Profiling a whole read set
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ProfileResult:
    panel: ReferencePanel
    reads: tuple                       # tuple[ReadCall, ...]

    def __post_init__(self) -> None:
        # Cache the completed subset once; it is read ~8x during reporting and the
        # dataclass is frozen, so set through object.__setattr__.
        object.__setattr__(
            self, "_complete",
            tuple(r for r in self.reads if r.status == "complete"),
        )

    @property
    def n_reads(self) -> int:
        return len(self.reads)

    @property
    def complete(self) -> tuple:
        return self._complete

    def status_counts(self) -> Counter:
        return Counter(r.status for r in self.reads)

    def composition(self) -> list:
        """[(genotype_tuple, n_reads)] over complete reads, most abundant first."""
        comp = Counter(r.genotype for r in self.complete)
        return comp.most_common()

    def domain_usage(self) -> dict:
        """domain_name -> Counter(source -> n_reads) over complete reads."""
        usage = {d.name: Counter() for d in self.panel.domains}
        for r in self.complete:
            for d in self.panel.domains:
                usage[d.name][r.calls[d.name].source] += 1
        return usage


def thin_marker_warnings(panel: ReferencePanel, min_markers: int) -> list[str]:
    """Sources whose private-marker pool in a domain is below the calling floor.

    Such a source can NEVER reach `min_markers` hits there, so every read carrying
    it in that domain is forced to 'partial' and silently dropped from composition.
    Surfacing this is load-bearing: otherwise the table reads 0% for a source that
    is actually present, with no hint why. Caller decides how to present it.
    """
    warnings: list[str] = []
    for d in panel.domains:
        for s in panel.sources:
            pool = panel.marker_pool[d.name][s]
            if pool < min_markers:
                warnings.append(
                    f"source '{s}' has only {pool} private {panel.k}-mers in domain "
                    f"'{d.name}' (< min_markers={min_markers}); reads using it there "
                    "will be counted 'partial', not in composition. Lower --min-markers "
                    "or --k, or check whether these sources are distinguishable here."
                )
    return warnings


def profile_reads(panel: ReferencePanel, records, *,
                  min_markers: int = DEFAULT_MIN_MARKERS,
                  margin: float = DEFAULT_MARGIN) -> ProfileResult:
    calls = tuple(
        call_read(rec.id, str(rec.seq), panel, min_markers=min_markers, margin=margin)
        for rec in records
    )
    return ProfileResult(panel, calls)


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def write_per_read_tsv(res: ProfileResult, path: Path) -> None:
    doms = [d.name for d in res.panel.domains]
    header = ["read_id", "length", "strand", "status"]
    for d in doms:
        header += [f"{d}_call", f"{d}_hits", f"{d}_frac"]
    lines = ["\t".join(header)]
    for r in res.reads:
        row = [r.read_id, str(r.length), r.strand, r.status]
        for d in doms:
            c = r.calls[d]
            row += [c.source or "-", str(c.hits), f"{c.fraction:.3f}"]
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_composition_tsv(res: ProfileResult, path: Path) -> None:
    doms = [d.name for d in res.panel.domains]
    total = len(res.complete)
    lines = ["\t".join(doms + ["n_reads", "fraction"])]
    for geno, n in res.composition():
        frac = n / total if total else 0.0
        lines.append("\t".join(list(geno) + [str(n), f"{frac:.4f}"]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_domain_usage_tsv(res: ProfileResult, path: Path) -> None:
    total = len(res.complete)
    usage = res.domain_usage()
    lines = ["\t".join(["domain", "source", "n_reads", "fraction"])]
    for d in res.panel.domains:
        for s in res.panel.sources:
            n = usage[d.name][s]
            frac = n / total if total else 0.0
            lines.append("\t".join([d.name, s, str(n), f"{frac:.4f}"]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_domain_usage(res: ProfileResult, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sources = res.panel.sources
    doms = [d.name for d in res.panel.domains]
    total = max(len(res.complete), 1)
    usage = res.domain_usage()
    grid = [[usage[d][s] / total for d in doms] for s in sources]

    fig, ax = plt.subplots(figsize=(1.6 + 1.1 * len(doms), 1.2 + 0.42 * len(sources)))
    im = ax.imshow(grid, aspect="auto", cmap="viridis", vmin=0)
    ax.set_xticks(range(len(doms)), doms)
    ax.set_yticks(range(len(sources)), sources)
    ax.set_xlabel("domain")
    ax.set_title(f"Source usage per domain  (n={len(res.complete)} complete reads)")
    for yi in range(len(sources)):
        for xi in range(len(doms)):
            v = grid[yi][xi]
            ax.text(xi, yi, f"{v*100:.0f}", ha="center", va="center",
                    color="white" if v < 0.6 * max(max(r) for r in grid) else "black",
                    fontsize=8)
    fig.colorbar(im, ax=ax, label="fraction of reads")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_top_genotypes(res: ProfileResult, path: Path, top: int = 20) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    full = res.composition()
    comp = full[:top]
    if not comp:
        return
    labels = [" | ".join(g) for g, _ in comp][::-1]
    counts = [n for _, n in comp][::-1]
    fig, ax = plt.subplots(figsize=(8, 1.0 + 0.32 * len(labels)))
    ax.barh(range(len(labels)), counts, color="#3b6ea5")
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.set_xlabel("reads")
    ax.set_title(f"Top {len(comp)} source-combinations "
                 f"(of {len(full)} observed)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _read_expected(path: Path, n_domains: int) -> set:
    """Optional designed-combination list: one combo per line, domains tab/comma
    separated, in the same order as detected domains."""
    combos: set = set()
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.replace(",", "\t").split("\t") if p.strip()]
        if len(parts) == n_domains:
            combos.add(tuple(parts))
    return combos


def build_summary(res: ProfileResult, ref_name: str, expected: set | None = None) -> str:
    doms = [d.name for d in res.panel.domains]
    sc = res.status_counts()
    total = res.n_reads
    n_complete = len(res.complete)
    comp = res.composition()
    usage = res.domain_usage()

    out: list[str] = []
    out.append(f"Chimera library profile  |  references: {ref_name} "
               f"({len(res.panel.sources)} sources, {len(doms)} domains: {', '.join(doms)})")
    out.append("=" * 72)
    out.append(f"reads processed         : {total}")
    out.append(f"  complete genotypes    : {n_complete} "
               f"({100*n_complete/total:.1f}%)" if total else "  complete genotypes    : 0")
    out.append(f"  partial (low coverage): {sc.get('partial', 0)}  "
               "<- truncated reads, or a domain source NOT in your panel")
    out.append(f"  ambiguous             : {sc.get('ambiguous', 0)}  "
               "<- no source cleared the margin (mixed/low-quality read)")
    out.append("")
    out.append(f"unique combinations observed : {len(comp)}")
    # With no designed list, the natural denominator is the panel's own
    # combinatorial space: (n_sources ^ n_domains) if every source can sit in
    # every domain. Real designs are usually a subset, so pass --expected for
    # the true coverage; this is an upper-bound frame.
    if expected is None and n_complete:
        space = len(res.panel.sources) ** len(doms)
        out.append(f"  of panel space        : {len(comp)} / {space} "
                   f"({100*len(comp)/space:.1f}%)  "
                   f"[{len(res.panel.sources)} sources ^ {len(doms)} domains; "
                   "upper bound -- use --expected for the designed subset]")

    # Per-domain dropouts: sources in the panel never used at that domain.
    out.append("")
    out.append("per-domain source usage (complete reads):")
    for d in res.panel.domains:
        used = usage[d.name]
        seen = [s for s in res.panel.sources if used[s] > 0]
        missing = [s for s in res.panel.sources if used[s] == 0]
        out.append(f"  {d.name}: {len(seen)}/{len(res.panel.sources)} sources seen"
                   + (f"; DROPOUT (never seen): {', '.join(missing)}" if missing else ""))

    # Abundance skew.
    if comp:
        top_g, top_n = comp[0]
        out.append("")
        out.append("abundance / skew:")
        out.append(f"  most abundant : {' | '.join(top_g)}  "
                   f"({top_n} reads, {100*top_n/n_complete:.1f}% of complete)")
        out.append(f"  fold-skew (top vs rarest observed) : "
                   f"{top_n / comp[-1][1]:.0f}x")

    # Optional coverage of the designed space.
    if expected is not None:
        observed = {g for g, _ in comp}
        hit = observed & expected
        out.append("")
        out.append("designed-space coverage:")
        out.append(f"  designed combinations : {len(expected)}")
        out.append(f"  observed of designed  : {len(hit)} "
                   f"({100*len(hit)/len(expected):.1f}%)")
        off = observed - expected
        out.append(f"  off-design observed   : {len(off)} "
                   "(combinations not in the designed list)")

    out.append("")
    out.append("CAVEATS (read before trusting counts):")
    out.append("  - Read counts are NOT absolute abundances: Nanopore yield is")
    out.append("    length- and GC-biased, and shorter molecules amplify/load better.")
    out.append("  - 'partial' reads are excluded from composition by design; a high")
    out.append("    partial rate can mean truncated reads OR sources missing from your")
    out.append("    reference panel. Add the missing source references to resolve them.")
    out.append("  - A combination being present is evidence of assembly, not of function.")
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _gather_reads(args: list[str]) -> list[Path]:
    paths: list[Path] = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            for ext, fmt in _READ_FORMATS.items():
                if fmt:
                    paths.extend(sorted(p.glob(f"*{ext}")))
        elif p.suffix.lower() in _READ_FORMATS and p.exists():
            paths.append(p)
    return paths


def _iter_records(paths: list[Path]):
    for p in paths:
        fmt = _READ_FORMATS.get(p.suffix.lower()) or "fastq"
        yield from SeqIO.parse(str(p), fmt)


class _Opts(TypedDict):
    k: int
    min_markers: int
    margin: float
    anchor_min: int
    names: tuple[str, ...] | None
    expected: Path | None
    positional: list[str]


def _nextarg(it, flag: str) -> str:
    """Pull the value following an option, with a friendly error if it's missing."""
    try:
        return next(it)
    except StopIteration:
        raise SystemExit(f"ERROR: option '{flag}' requires a value.")


def _parse_args(argv: list[str]) -> _Opts:
    opts: _Opts = {"k": DEFAULT_K, "min_markers": DEFAULT_MIN_MARKERS,
                   "margin": DEFAULT_MARGIN, "anchor_min": DEFAULT_ANCHOR_MIN,
                   "names": None, "expected": None, "positional": []}
    it = iter(argv)
    for tok in it:
        if tok == "--k":
            opts["k"] = int(_nextarg(it, tok))
        elif tok == "--min-markers":
            opts["min_markers"] = int(_nextarg(it, tok))
        elif tok == "--margin":
            opts["margin"] = float(_nextarg(it, tok))
        elif tok == "--anchor-min":
            opts["anchor_min"] = int(_nextarg(it, tok))
        elif tok == "--names":
            opts["names"] = tuple(s.strip() for s in _nextarg(it, tok).split(","))
        elif tok == "--expected":
            opts["expected"] = Path(_nextarg(it, tok))
        elif tok.startswith("--"):
            raise SystemExit(f"ERROR: unknown option '{tok}'")
        else:
            opts["positional"].append(tok)
    return opts


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 1
    opts = _parse_args(argv)
    pos = opts["positional"]
    if len(pos) < 2:
        print("ERROR: need a references FASTA and at least one reads file/folder.")
        return 1

    ref_path = Path(pos[0])
    if not ref_path.exists():
        print(f"ERROR: references not found: {ref_path}")
        return 1
    try:
        refs = read_references(ref_path)
    except Exception as e:
        print(f"ERROR: could not read references '{ref_path.name}': {e}")
        return 1

    read_paths = _gather_reads(pos[1:])
    if not read_paths:
        print("ERROR: no reads found. Pass .fastq/.fasta files or a folder.")
        return 1

    print(f"Building panel from {len(refs)} references ({ref_path.name}) ...")
    try:
        panel = build_panel(refs, k=opts["k"], anchor_min=opts["anchor_min"],
                            names=opts["names"])
    except Exception as e:
        print(f"ERROR: {e}")
        return 1
    print(f"  detected {len(panel.domains)} domains: "
          + ", ".join(f"{d.name}[{d.start}-{d.end}]" for d in panel.domains))
    print("  private markers / domain: "
          + ", ".join(f"{d.name}:{sum(panel.marker_pool[d.name].values())}"
                      for d in panel.domains))
    for w in thin_marker_warnings(panel, opts["min_markers"]):
        print(f"  WARNING: {w}")

    records = _iter_records(read_paths)
    res = profile_reads(panel, records, min_markers=opts["min_markers"],
                        margin=opts["margin"])
    print(f"Profiled {res.n_reads} reads from {len(read_paths)} file(s).")

    expected = None
    if opts["expected"] is not None:
        if opts["expected"].exists():
            expected = _read_expected(opts["expected"], len(panel.domains))
        else:
            print(f"  NOTE: --expected file not found: {opts['expected']}")

    out_dir = ref_path.parent / "library_profile"
    out_dir.mkdir(exist_ok=True)
    write_per_read_tsv(res, out_dir / "per_read.tsv")
    write_composition_tsv(res, out_dir / "composition.tsv")
    write_domain_usage_tsv(res, out_dir / "domain_usage.tsv")
    plot_domain_usage(res, out_dir / "domain_usage.png")
    plot_top_genotypes(res, out_dir / "top_genotypes.png")
    summary = build_summary(res, ref_path.stem, expected)
    (out_dir / "library_summary.txt").write_text(summary, encoding="utf-8")

    print()
    print(summary)
    print(f"Outputs -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
