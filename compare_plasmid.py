"""Compare a whole-plasmid sequence (e.g. Plasmidsaurus) to a SnapGene .dna map.

Companion to align_sanger.py. Where that script aligns short Sanger .ab1 reads
to a region of a map, THIS one compares a full assembled plasmid (the kind you
get back from Plasmidsaurus Nanopore whole-plasmid sequencing) against your
reference map and tells you whether they are identical, and if not, exactly
where and how they differ.

It handles the two things that make whole-plasmid comparison different from a
Sanger read:
  - CIRCULARITY: an assembly can start at any base around the circle, so the
    reference is "doubled" internally and differences are reported in your
    map's coordinates regardless of where the assembly happens to begin.
  - STRAND: the assembly may be the reverse complement of your map; both
    orientations are tried automatically.

There is no per-base quality in a consensus FASTA, so (unlike the Sanger
script) every difference is reported - a polished consensus is already
high-confidence.

HOW IT WORKS (the pipeline, in order)
    1. read the .dna map (reference) and the .fasta/.gbk whole-plasmid (query).
    2. DOUBLE the reference (ref + ref). A circular assembly can be "cut" at
       any base, so any rotation of it appears as one continuous stretch
       somewhere inside the doubled sequence -> see compare(), line ~161.
    3. Local-align the query AND its reverse complement to the doubled
       reference; keep whichever strand scores higher -> compare(), ~165-172.
       Gaps are penalised heavily so real point mutations are called as
       mismatches, not hidden behind invented gaps -> _make_aligner().
    4. WALK the alignment column by column, classifying each as match,
       mismatch, insertion, or deletion. Positions are taken modulo the real
       reference length, so every difference is reported in YOUR map's
       1-based coordinates no matter where the assembly started -> _walk().
    5. ANNOTATE each difference using the map's features: is it in a coding
       sequence (and is it synonymous / missense / nonsense), or in a backbone
       element like an origin, promoter, terminator, or resistance gene?
       -> annotate_difference(); a one-line verdict comes from _assess().
    6. Report it three ways: console text, a figure, and a summary file.
    Tip: every function is importable and runs standalone, so you can call
    read_reference / read_query / compare in a Python prompt to inspect a
    result while troubleshooting.

USAGE
    python compare_plasmid.py REFERENCE.dna QUERY.fasta [MORE.fasta/.gbk ...]
    python compare_plasmid.py REFERENCE.dna  FOLDER_OF_RESULTS

OUTPUTS (in "comparison_figures/" next to the reference .dna)
    <query>__vs__<ref>.png      one figure per query
    comparison_summary.txt       one-line IDENTICAL / differences verdict per query

Dependencies: biopython, matplotlib   (see README.txt)
"""
from __future__ import annotations

import sys
import textwrap
from dataclasses import dataclass, replace
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.Seq import Seq
except ModuleNotFoundError as _exc:   # friendly first-run message for non-programmers
    raise SystemExit(
        f"ERROR: missing required package '{_exc.name}'.\n"
        f"  Install dependencies once with:  pip install -r requirements.txt\n"
        f"  (or:  pip install biopython matplotlib)"
    )

# --- Scoring: tuned for high-identity whole-plasmid comparison. ------------
MATCH, MISMATCH = 2.0, -3.0
OPEN_GAP, EXTEND_GAP = -7.0, -2.0

CONTEXT_BP = 12               # bases of sequence context shown around a difference
LARGE_PLASMID_WARN = 30_000   # warn above this size: alignment gets slow (O(N^2) time).
                              # NOT a memory limit - Biopython 1.87 aligns in linear
                              # space, so RSS stays flat (measured, AUDIT_REPORT ISSUE-2).
TESTED_BIOPYTHON = "1.87"     # version the suite + memory-safety were validated against

# Map common extensions to the Biopython format name.
# NOTE on .ab1: Plasmidsaurus sometimes delivers the whole-plasmid CONSENSUS as
# an .ab1 (so it opens in trace viewers). Here we treat an .ab1 query as a full
# plasmid sequence (quality ignored). That is different from align_sanger.py,
# which treats an .ab1 as a short Sanger READ covering only a region of the map.
# Rule of thumb: full-length .ab1 (~= plasmid size) -> this script; short
# (~1 kb) Sanger read .ab1 -> align_sanger.py.
_QUERY_FORMATS = {
    ".fasta": "fasta", ".fa": "fasta", ".fan": "fasta", ".fna": "fasta", ".seq": "fasta",
    ".gb": "genbank", ".gbk": "genbank", ".genbank": "genbank",
    ".ab1": "abi",
}


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def warn_if_untested_biopython() -> None:
    """Warn (do not block) when Biopython differs from the validated version.

    CONCERN-5: load-bearing behaviors are version-dependent - linear-space alignment
    (memory-safety), translation, minus-strand FeatureLocation iteration, and
    aligner coordinates. A different version may silently change any of them.
    """
    import Bio
    if Bio.__version__ != TESTED_BIOPYTHON:
        print(f"NOTE: this toolkit was validated against Biopython {TESTED_BIOPYTHON}; "
              f"you have {Bio.__version__}. Re-run test_toolkit.py to confirm "
              f"alignment/translation behavior (and memory-safety) are unchanged.")


def read_reference(dna_path: Path):
    """Read a reference map. Accepts a SnapGene .dna or a GenBank .gb/.gbk file
    (both carry the CDS/feature annotations the analysis needs)."""
    fmt = "genbank" if dna_path.suffix.lower() in (".gb", ".gbk", ".gbff", ".genbank") else "snapgene"
    rec = SeqIO.read(str(dna_path), fmt)
    if not rec.name or rec.name == "<unknown name>":
        rec.name = dna_path.stem
    return rec


def read_query(path: Path):
    """Read a whole-plasmid sequence from FASTA or GenBank.

    A consensus is normally a single record. If a file holds several contigs
    (an incomplete or multimeric assembly), use the longest and note the rest
    rather than skipping the whole file (which SeqIO.read would force).
    """
    fmt = _QUERY_FORMATS.get(path.suffix.lower())
    if fmt is None:
        raise ValueError(f"unsupported query format '{path.suffix}' (use .fasta or .gbk)")
    records = list(SeqIO.parse(str(path), fmt))
    if not records:
        raise ValueError(f"no sequences found in {path.name}")
    if len(records) > 1:
        records.sort(key=lambda r: len(r.seq), reverse=True)
        dropped = ", ".join(f"{r.id} ({len(r.seq):,} bp)" for r in records[1:])
        print(f"  NOTE: {path.name} has {len(records)} contigs; using the longest "
              f"({len(records[0].seq):,} bp), ignoring: {dropped}")
    rec = records[0]
    # Always label outputs by the FILENAME, not the record's internal id.
    # (A Plasmidsaurus FASTA header can be a short tag like "5.1-1" that does
    # not match the file name - using the file name avoids confusion.)
    rec.id = path.stem
    return rec


# --------------------------------------------------------------------------- #
# Comparison
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Difference:
    ref_pos: int      # 1-based position in the reference map (start, for indels)
    kind: str         # "mismatch" | "insertion" | "deletion"
    ref_base: str     # may be multi-base after consecutive deletions are collapsed
    query_base: str   # may be multi-base after consecutive insertions are collapsed
    col: int          # column in the gapped alignment (for context display)

    @property
    def span_len(self) -> int:
        """Number of bases involved (deleted, inserted, or 1 for a mismatch)."""
        if self.kind == "deletion":
            return len(self.ref_base)
        if self.kind == "insertion":
            return len(self.query_base)
        return 1

    @property
    def ref_range(self) -> str:
        """Human 1-based reference span, e.g. '4,000' or '4,000-4,002'."""
        if self.kind == "deletion" and self.span_len > 1:
            return f"{self.ref_pos:,}-{self.ref_pos + self.span_len - 1:,}"
        return f"{self.ref_pos:,}"


@dataclass(frozen=True)
class Consequence:
    """Biological meaning of a difference, derived from the map's annotations."""
    category: str     # synonymous | missense | nonsense | stop-lost |
                      # coding-indel | feature | intergenic
    text: str         # human-readable, e.g. "missense p.Thr46Ala in CDS 'MyGene'"
    feature: str      # feature label it falls in (or '')

    @property
    def is_coding_change(self) -> bool:
        return self.category in {"synonymous", "missense", "nonsense",
                                 "stop-lost", "coding-indel"}

    @property
    def is_protein_altering(self) -> bool:
        return self.category in {"missense", "nonsense", "stop-lost", "coding-indel"}


@dataclass(frozen=True)
class CompareResult:
    query_name: str
    ref_name: str
    strand: str               # "+" or "-"
    ref_len: int
    query_len: int
    aligned_length: int
    matches: int
    differences: tuple[Difference, ...]
    consequences: tuple[Consequence, ...]   # parallel to differences
    assessment: str                         # one-line biological summary
    aligned_ref: str          # gapped reference row
    aligned_query: str        # gapped query row
    n_lowconf: int            # lowercase ("low-confidence") bases seen in the query
    clipped_query_bp: int     # query bases soft-clipped at the ends (not checked)
    examined_ref_bp: int      # reference bases the alignment actually traversed

    @property
    def percent_identity(self) -> float:
        return 100.0 * self.matches / self.aligned_length if self.aligned_length else 0.0

    @property
    def n_mismatch(self) -> int:
        return sum(1 for d in self.differences if d.kind == "mismatch")

    @property
    def n_indel(self) -> int:
        return sum(1 for d in self.differences if d.kind != "mismatch")

    @property
    def is_identical(self) -> bool:
        return not self.differences

    @property
    def is_fully_examined(self) -> bool:
        """True only when a clean verdict would be EARNED: the whole query was used
        (no soft-clip), the entire reference circle was covered, and no low-confidence
        bases were present. If any is false, a green all-clear is not justified."""
        return (self.clipped_query_bp == 0
                and self.examined_ref_bp >= self.ref_len
                and self.n_lowconf == 0)

    @property
    def completeness_caveats(self) -> tuple[str, ...]:
        """Why the molecule was not fully examined (empty when it was)."""
        caveats = []
        if self.clipped_query_bp:
            caveats.append(f"{self.clipped_query_bp} bp of query soft-clipped/unexamined")
        if self.examined_ref_bp < self.ref_len:
            caveats.append(f"{self.ref_len - self.examined_ref_bp} bp of the map not covered")
        if self.n_lowconf:
            caveats.append(f"{self.n_lowconf} low-confidence base(s)")
        return tuple(caveats)

    @property
    def verdict(self) -> str:
        """Honest one-line verdict, gated on verified completeness."""
        if not self.is_identical:
            return f"{len(self.differences)} difference(s): {self.n_mismatch} mismatch, {self.n_indel} indel"
        if self.is_fully_examined:
            return "IDENTICAL to map"
        return ("NO DIFFERENCES in the examined region - NOT a verified all-clear ("
                + "; ".join(self.completeness_caveats) + ")")


def _make_aligner() -> PairwiseAligner:
    aln = PairwiseAligner()
    aln.mode = "local"
    aln.match_score = MATCH
    aln.mismatch_score = MISMATCH
    aln.open_gap_score = OPEN_GAP
    aln.extend_gap_score = EXTEND_GAP
    return aln


def _walk(alignment, ref_offset: int, ref_len: int):
    """Walk gapped columns; report differences in 1-based reference coordinates.

    ref_offset is the start in the DOUBLED reference; we take positions modulo
    the real reference length so coordinates land in the original map's frame.
    """
    target = alignment[0]   # (doubled) reference row, '-' for gaps
    query = alignment[1]    # query row, '-' for gaps
    ref_idx = ref_offset
    diffs: list[Difference] = []
    matches = 0
    for col, (t, q) in enumerate(zip(target, query)):
        if t == "-":                          # extra base in query (insertion)
            diffs.append(Difference(ref_idx % ref_len + 1, "insertion", "-", q, col))
        elif q == "-":                        # base missing from query (deletion)
            diffs.append(Difference(ref_idx % ref_len + 1, "deletion", t, "-", col))
            ref_idx += 1
        else:
            if t.upper() == q.upper():
                matches += 1
            else:
                diffs.append(Difference(ref_idx % ref_len + 1, "mismatch", t, q, col))
            ref_idx += 1
    return tuple(diffs), matches


def _collapse_indels(diffs: tuple[Difference, ...]) -> tuple[Difference, ...]:
    """Merge runs of consecutive single-base indels into one multi-base event.

    A 3 bp deletion comes out of _walk() as three adjacent single-base deletion
    columns; this joins them into one "3 bp deletion". Adjacency is judged by
    alignment column, so only a truly contiguous gap is merged. Mismatches are
    never merged. Differences arrive in column order, so a single pass suffices.
    """
    out: list[Difference] = []
    for d in diffs:
        if out and d.kind == out[-1].kind and d.kind in ("insertion", "deletion"):
            prev = out[-1]
            prev_cols = len(prev.ref_base) if d.kind == "deletion" else len(prev.query_base)
            if d.col == prev.col + prev_cols:          # exactly the next column
                if d.kind == "deletion":
                    out[-1] = replace(prev, ref_base=prev.ref_base + d.ref_base)
                else:
                    out[-1] = replace(prev, query_base=prev.query_base + d.query_base)
                continue
        out.append(d)
    return tuple(out)


# --------------------------------------------------------------------------- #
# Biological annotation - what does each difference MEAN?
# --------------------------------------------------------------------------- #
_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu",
    "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "*", "X": "Xaa",
}
# Feature types that don't help locate a difference (whole-plasmid spans).
_IGNORED_FEATURE_TYPES = {"source"}


def _flabel(f) -> str:
    """Best human-readable label for a feature."""
    for key in ("label", "gene", "product", "note"):
        if key in f.qualifiers:
            return f.qualifiers[key][0]
    return f.type


def _overlapping(ref, pos0: int):
    """Annotated features covering a 0-based reference position."""
    return [f for f in ref.features
            if f.type not in _IGNORED_FEATURE_TYPES and pos0 in f.location]


def _codon_consequence(ref, cds, pos0: int, query_base_fwd: str) -> Consequence | None:
    """For a mismatch inside a CDS, translate the affected codon both ways.

    query_base_fwd is the query base on the FORWARD reference strand (that is
    how the alignment stores it). For a minus-strand CDS we complement it so
    the substitution is applied in coding orientation.
    """
    if len(query_base_fwd) != 1 or query_base_fwd in "-Nn":
        return None
    # extract() does NOT apply a GenBank codon_start offset, so a CDS starting on
    # a partial codon (codon_start != 1) would frame-shift every residue call.
    # SnapGene full-ORF maps always use codon_start=1; guard rather than mis-call.
    if int(cds.qualifiers.get("codon_start", [1])[0]) != 1:
        lbl = _flabel(cds)
        return Consequence("feature",
                           f"in CDS '{lbl}' (codon_start != 1; protein effect not called)", lbl)
    coding_positions = list(cds.location)          # genomic positions, coding order
    try:
        ci = coding_positions.index(pos0)          # index into the coding sequence
    except ValueError:
        return None
    coding = str(cds.extract(ref.seq)).upper()     # coding-sense DNA (strand-aware)
    codon_start = (ci // 3) * 3
    if codon_start + 3 > len(coding):
        return None
    ref_codon = coding[codon_start:codon_start + 3]
    new_base = query_base_fwd.upper()
    if cds.location.strand == -1:
        new_base = str(Seq(new_base).complement())
    within = ci % 3
    new_codon = ref_codon[:within] + new_base + ref_codon[within + 1:]
    try:
        aa_ref = str(Seq(ref_codon).translate())
        aa_new = str(Seq(new_codon).translate())
    except Exception:
        return None

    aa_num = codon_start // 3 + 1                   # 1-based residue number
    lbl = _flabel(cds)
    r3, n3 = _AA3.get(aa_ref, aa_ref), _AA3.get(aa_new, aa_new)
    if aa_ref == aa_new:
        return Consequence("synonymous", f"synonymous p.{r3}{aa_num}= in CDS '{lbl}'", lbl)
    if aa_new == "*":
        return Consequence("nonsense", f"NONSENSE p.{r3}{aa_num}* in CDS '{lbl}' (premature stop)", lbl)
    if aa_ref == "*":
        return Consequence("stop-lost", f"stop-lost p.*{aa_num}{n3} in CDS '{lbl}'", lbl)
    return Consequence("missense", f"missense p.{r3}{aa_num}{n3} in CDS '{lbl}'", lbl)


def annotate_difference(ref, d: Difference) -> Consequence:
    """Classify one difference by what map feature it hits and its protein effect."""
    pos0 = d.ref_pos - 1
    overlapping = _overlapping(ref, pos0)
    cds_list = [f for f in overlapping if f.type == "CDS"]

    if d.kind == "mismatch":
        for cds in cds_list:
            con = _codon_consequence(ref, cds, pos0, d.query_base)
            if con:
                return con
    elif cds_list:                                 # insertion / deletion inside a CDS
        lbl = _flabel(cds_list[0])
        n = d.span_len
        frame = "in-frame (multiple of 3)" if n % 3 == 0 else "FRAMESHIFT"
        return Consequence("coding-indel",
                           f"{n} bp {d.kind} in CDS '{lbl}' - {frame}", lbl)

    if overlapping:                                # non-coding feature (ori, promoter, ...)
        f = min(overlapping, key=lambda x: int(x.location.end) - int(x.location.start))
        lbl = _flabel(f)
        return Consequence("feature", f"in {f.type} '{lbl}'", lbl)
    return Consequence("intergenic", "intergenic (no annotated feature)", "")


def _assess(consequences: tuple[Consequence, ...], fully_examined: bool) -> str:
    """One-line biological verdict across all differences for a query.

    A green all-clear ("IDENTICAL", "ORFs intact") is only emitted when the
    molecule was fully examined; otherwise the claim is qualified to the examined
    region so the verdict never over-promises (see CONCERN-1 / is_fully_examined).
    """
    if not consequences:
        if fully_examined:
            return "IDENTICAL to map - no differences."
        return ("NO DIFFERENCES in the examined region - NOT a verified all-clear "
                "(unexamined region or low-confidence bases; see caveats).")
    altering = [c for c in consequences if c.is_protein_altering]  # missense/nonsense/etc.
    if altering:
        return "ATTENTION - protein-altering change(s): " + "; ".join(c.text for c in altering)
    # Only synonymous and/or non-coding differences remain.
    locs = sorted({c.feature or "intergenic" for c in consequences})
    silent = sum(1 for c in consequences if c.category == "synonymous")
    note = f" ({silent} synonymous)" if silent else ""
    intact = "insert/ORFs intact" if fully_examined else "insert/ORFs intact IN THE EXAMINED REGION (not fully verified)"
    return f"No protein-coding changes - {intact}{note}. Difference(s) in: {', '.join(locs)}."


def compare(ref, query) -> CompareResult:
    """Align a whole-plasmid query to a circular reference; report differences."""
    aligner = _make_aligner()
    ref_len = len(ref.seq)
    # Double the reference so an assembly starting anywhere on the circle (or
    # spanning the origin) still aligns as one contiguous block.
    doubled = (str(ref.seq) + str(ref.seq)).upper()
    n_lowconf = sum(1 for b in str(query.seq) if b.islower())
    q_up = str(query.seq).upper()

    candidates = []
    for strand, seq in (("+", q_up), ("-", str(query.seq.reverse_complement()).upper())):
        if seq:
            candidates.append((strand, aligner.align(doubled, seq)[0]))
    if not candidates:
        raise ValueError("query has no sequence")

    strand, aln = max(candidates, key=lambda c: c[1].score)
    ref_offset = int(aln.coordinates[0][0])
    # Local alignment can soft-clip query bases at the termini; a real difference
    # within ~1 bp of the (arbitrary) cut point of the circular consensus would
    # then be dropped and never lower identity. Count clipped query bases so the
    # report can warn. (A large count usually means the wrong reference map.)
    q_start, q_end = int(aln.coordinates[1][0]), int(aln.coordinates[1][-1])
    clipped_query_bp = q_start + (len(query.seq) - q_end)
    # Reference bases the alignment actually traversed; < ref_len means part of the
    # map (and any CDS there) was never examined, so a clean verdict isn't earned.
    examined_ref_bp = int(aln.coordinates[0][-1]) - int(aln.coordinates[0][0])
    raw_differences, matches = _walk(aln, ref_offset, ref_len)
    # % identity is a per-BASE measure, so compute the length before collapsing.
    aligned_length = matches + len(raw_differences)
    differences = _collapse_indels(raw_differences)
    consequences = tuple(annotate_difference(ref, d) for d in differences)
    fully_examined = (clipped_query_bp == 0 and examined_ref_bp >= ref_len and n_lowconf == 0)

    return CompareResult(
        query_name=query.id,
        ref_name=ref.name,
        strand=strand,
        ref_len=ref_len,
        query_len=len(query.seq),
        aligned_length=aligned_length,
        matches=matches,
        differences=differences,
        consequences=consequences,
        assessment=_assess(consequences, fully_examined),
        aligned_ref=str(aln[0]),
        aligned_query=str(aln[1]),
        n_lowconf=n_lowconf,
        clipped_query_bp=clipped_query_bp,
        examined_ref_bp=examined_ref_bp,
    )


# --------------------------------------------------------------------------- #
# Sequence-context display
# --------------------------------------------------------------------------- #
def context_lines(res: CompareResult, d: Difference) -> tuple[str, str, str]:
    """Return (ref, match-bar, query) text around a difference, gaps aligned."""
    lo = max(0, d.col - CONTEXT_BP)
    hi = d.col + CONTEXT_BP + 1
    ref_chunk = res.aligned_ref[lo:hi]
    qry_chunk = res.aligned_query[lo:hi]
    bar = "".join(
        "|" if a != "-" and b != "-" and a.upper() == b.upper() else " "
        for a, b in zip(ref_chunk, qry_chunk)
    )
    return ref_chunk, bar, qry_chunk


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def print_report(res: CompareResult) -> None:
    print(f"\n{'='*70}")
    print(f"  {res.query_name}")
    print(f"  vs {res.ref_name}  (strand {res.strand})")
    print(f"{'='*70}")
    print(f"  Identity : {res.percent_identity:.2f}%  ({res.matches}/{res.aligned_length})")
    print(f"  Length   : query {res.query_len:,} bp  vs  map {res.ref_len:,} bp")
    print(f"  Verdict  : {res.verdict}")
    print(f"  Assess   : {res.assessment}")
    if res.n_lowconf:
        print(f"  Note     : query has {res.n_lowconf} lowercase (low-confidence) base(s)")
    if res.clipped_query_bp:
        print(f"  WARNING  : {res.clipped_query_bp} query base(s) soft-clipped at the ends and NOT "
              f"checked -")
        print(f"             a difference there would be missed. A large count usually means the "
              f"wrong reference map.")
    for d, c in list(zip(res.differences, res.consequences))[:50]:
        label = {"mismatch": f"ref {d.ref_base} -> query {d.query_base}",
                 "insertion": f"{d.span_len} bp inserted in query: {d.query_base}",
                 "deletion": f"{d.span_len} bp missing from query (ref {d.ref_base})"}[d.kind]
        ref_c, bar, qry_c = context_lines(res, d)
        print(f"  {'-'*50}")
        print(f"  ref {d.ref_range}  {d.kind}: {label}")
        print(f"      -> {c.text}")
        print(f"      ref   5'-{ref_c}-3'")
        print(f"              {bar}")
        print(f"      query 5'-{qry_c}-3'")
    if len(res.differences) > 50:
        print(f"  ... and {len(res.differences) - 50} more")


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
_COLORS = {"mismatch": "#d62728", "insertion": "#1f77b4", "deletion": "#9467bd"}


def plot_comparison(res: CompareResult, out_path: Path) -> None:
    """Figure: full circular reference as a bar with difference ticks + context."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig, (ax_map, ax_txt) = plt.subplots(
        2, 1, figsize=(12, 5.5), height_ratios=[1, 1.8], gridspec_kw={"hspace": 0.4},
    )

    # Panel 1 - whole reference, with a tick at each difference.
    ax_map.add_patch(plt.Rectangle((1, 0.3), res.ref_len, 0.4, color="#2ca02c", alpha=0.2))
    ax_map.hlines(0.5, 1, res.ref_len, color="#2ca02c", lw=2)
    for d in res.differences:
        ax_map.vlines(d.ref_pos, 0.1, 0.9, color=_COLORS[d.kind], lw=1.3)
        if d.kind == "mismatch":
            tag = f"{d.ref_base}>{d.query_base}"
        else:
            tag = f"{d.kind[:3]}{d.span_len}" if d.span_len > 1 else d.kind[:3]
        ax_map.annotate(tag, (d.ref_pos, 0.93), fontsize=7, ha="center", color=_COLORS[d.kind])
    ax_map.set_xlim(-res.ref_len * 0.02, res.ref_len * 1.02)
    ax_map.set_ylim(0, 1.05)
    ax_map.set_yticks([])
    ax_map.set_xlabel("reference map position (bp)  -  circular")
    seen = {d.kind for d in res.differences}
    if seen:
        ax_map.legend(handles=[Patch(color=_COLORS[k], label=k) for k in sorted(seen)],
                      loc="upper right", fontsize=8, ncol=3)

    # Panel 2 - biological assessment + aligned sequence around each difference.
    ax_txt.axis("off")
    if res.is_identical and res.is_fully_examined:
        ax_txt.text(0.5, 0.6, "IDENTICAL to reference map", ha="center", va="center",
                    fontsize=15, color="#2ca02c", weight="bold")
        ax_txt.text(0.5, 0.3, f"{res.percent_identity:.2f}% over {res.aligned_length:,} bp",
                    ha="center", va="center", fontsize=10, color="#555")
    elif res.is_identical:
        # No differences, but the molecule was NOT fully examined - amber, not green.
        ax_txt.text(0.5, 0.62, "NO DIFFERENCES in the examined region", ha="center", va="center",
                    fontsize=13, color="#d9820a", weight="bold")
        ax_txt.text(0.5, 0.42, "NOT a verified all-clear", ha="center", va="center",
                    fontsize=11, color="#d9820a")
        ax_txt.text(0.5, 0.2, "; ".join(res.completeness_caveats), ha="center", va="center",
                    fontsize=8.5, color="#555")
    else:
        # Headline assessment, coloured by severity, wrapped to fit.
        protein_altering = any(c.is_protein_altering for c in res.consequences)
        head_color = "#d62728" if protein_altering else "#2ca02c"
        wrapped = textwrap.fill(res.assessment, width=95)
        ax_txt.text(0.01, 0.99, wrapped, ha="left", va="top", fontsize=9.5,
                    weight="bold", color=head_color, transform=ax_txt.transAxes)

        n_head = wrapped.count("\n") + 1
        lines = []
        for d, c in list(zip(res.differences, res.consequences))[:5]:
            ref_c, bar, qry_c = context_lines(res, d)
            head = f"ref {d.ref_range}  {d.kind}  ->  {c.text}"
            lines += [head, f"  ref   {ref_c}", f"        {bar}", f"  query {qry_c}", ""]
        if len(res.differences) > 5:
            lines.append(f"... and {len(res.differences) - 5} more (see console / summary)")
        # Place the monospace detail below the (1-2 line) headline.
        y_start = 0.99 - 0.075 * (n_head + 0.5)
        ax_txt.text(0.01, y_start, "\n".join(lines), ha="left", va="top",
                    fontfamily="monospace", fontsize=9, transform=ax_txt.transAxes)

    if res.is_identical:
        verdict = "IDENTICAL" if res.is_fully_examined else "NO DIFFS (not fully examined)"
    else:
        verdict = f"{res.n_mismatch} mismatch, {res.n_indel} indel"
    fig.suptitle(
        f"{res.query_name}  vs  {res.ref_name}  ({res.strand})\n"
        f"{res.percent_identity:.2f}% identity   |   query {res.query_len:,} bp / map {res.ref_len:,} bp"
        f"   |   {verdict}",
        fontsize=10, y=0.99,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _gather_queries(args: list[str]) -> list[Path]:
    queries: list[Path] = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            for ext in _QUERY_FORMATS:
                queries.extend(sorted(p.glob(f"*{ext}")))
        elif p.suffix.lower() in _QUERY_FORMATS and p.exists():
            queries.append(p)
    return queries


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 1
    warn_if_untested_biopython()

    ref_path = Path(argv[0])
    if not ref_path.exists():
        print(f"ERROR: reference not found: {ref_path}")
        return 1
    try:
        ref = read_reference(ref_path)
    except Exception as e:
        print(f"ERROR: could not read reference '{ref_path.name}' as a SnapGene .dna file: {e}")
        return 1

    queries = _gather_queries(argv[1:])
    if not queries:
        print("ERROR: no query sequences found. Pass .fasta/.gbk files or a folder.")
        return 1

    if len(ref.seq) > LARGE_PLASMID_WARN:
        print(f"NOTE: reference is large ({len(ref.seq):,} bp); alignment may be slow.")

    out_dir = ref_path.parent / "comparison_figures"
    out_dir.mkdir(exist_ok=True)
    print(f"Reference: {ref.name} ({len(ref.seq):,} bp, {ref.annotations.get('topology', '?')})")
    print(f"Comparing {len(queries)} quer(ies); outputs -> {out_dir}")

    summary: list[str] = []
    for qpath in queries:
        try:
            query = read_query(qpath)
        except Exception as e:           # one bad file must not stop the batch
            print(f"\n  SKIPPED {qpath.name}: {e}")
            summary.append(f"SKIPPED   {qpath.name}: {e}")
            continue
        if len(query.seq) < 0.5 * len(ref.seq):
            print(f"\n  WARNING: {qpath.name} is {len(query.seq):,} bp, much shorter than the "
                  f"{len(ref.seq):,} bp map.")
            print(f"           This script expects a WHOLE plasmid. If this is a short Sanger")
            print(f"           READ (~1 kb), use align_sanger.py instead.")
        res = compare(ref, query)
        print_report(res)
        # res.query_name is already the clean filename stem; do NOT call .stem
        # again (a name like "sample_5.1-1" would lose its ".1-1" tail).
        plot_comparison(res, out_dir / f"{_safe(res.query_name)}__vs__{_safe(ref.name)}.png")
        if res.is_identical:
            tag = "IDENTICAL" if res.is_fully_examined else "NO DIFFS (unverified)"
        else:
            tag = f"{res.n_mismatch} mismatch, {res.n_indel} indel"
        summary.append(f"{tag:<22} {res.percent_identity:6.2f}%  query {res.query_len:,} bp  {res.query_name}")
        summary.append(f"    -> {res.assessment}")
        for d, c in zip(res.differences, res.consequences):
            summary.append(f"       ref {d.ref_range} {d.kind}: {c.text}")

    summary_path = out_dir / "comparison_summary.txt"
    header = f"Whole-plasmid comparison  |  reference: {ref.name} ({len(ref.seq):,} bp)\n" + "=" * 70
    summary_path.write_text(header + "\n" + "\n".join(summary) + "\n", encoding="utf-8")
    print(f"\nSummary -> {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
