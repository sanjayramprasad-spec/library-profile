"""Align Sanger .ab1 sequencing reads against SnapGene .dna plasmid maps.

A simple replacement for the SnapGene "align to reference" workflow used for
cloning validation. It parses the proprietary .dna map and the .ab1 trace,
quality-trims the read, does a proper local (Smith-Waterman) alignment on both
strands, and reports % identity plus the number / location / sequence-context
of every real mismatch and indel, with figures.

USAGE
    python align_sanger.py REFERENCE.dna READ.ab1 [MORE.ab1 ...]
    python align_sanger.py REFERENCE.dna  FOLDER_OF_AB1S

OUTPUTS (written next to the reference, in a folder "alignment_figures/")
    <read>__vs__<ref>.png   one detail figure per read
    COMBINED__<ref>.png      all reads tiled on the reference (if >1 read)
    alignment_summary.txt    one-line PASS / discrepancy summary per read

Dependencies: biopython, matplotlib   (see README.txt)
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, replace
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
except ModuleNotFoundError as _exc:   # friendly first-run message for non-programmers
    raise SystemExit(
        f"ERROR: missing required package '{_exc.name}'.\n"
        f"  Install dependencies once with:  pip install -r requirements.txt\n"
        f"  (or:  pip install biopython matplotlib)"
    )

# --- Scoring: tuned for high-identity Sanger-vs-reference comparison. -------
# Gap penalties are stiff so we don't paper over real indels with cheap gaps.
MATCH, MISMATCH = 2.0, -3.0
OPEN_GAP, EXTEND_GAP = -7.0, -2.0

QUAL_TRIM_THRESHOLD = 20  # Phred; trim ends below this (Q20 = 99% base accuracy)
TRIM_WINDOW = 10          # sliding window (bp) for end trimming
MIN_CONF_Q = 30           # a discrepancy only "counts" if the read is this sure
                          # (Q30 = 99.9%); N / low-Q calls are unread, not mutations
CONTEXT_BP = 10           # bases of sequence context shown around each discrepancy
SANGER_MAX_BP = 1500      # a real Sanger read is < ~1.2 kb; longer means the file
                          # is probably a WHOLE-PLASMID result (use compare_plasmid.py)
TESTED_BIOPYTHON = "1.87"  # version the regression suite was validated against


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #
def warn_if_untested_biopython() -> None:
    """Warn (do not block) when Biopython differs from the validated version;
    translation, FeatureLocation, and aligner coordinates are version-dependent
    (CONCERN-5)."""
    import Bio
    if Bio.__version__ != TESTED_BIOPYTHON:
        print(f"NOTE: this toolkit was validated against Biopython {TESTED_BIOPYTHON}; "
              f"you have {Bio.__version__}. Re-run test_toolkit.py to confirm behavior.")


def read_reference(dna_path: Path):
    """Read a SnapGene .dna map. Returns a SeqRecord (with features, topology)."""
    rec = SeqIO.read(str(dna_path), "snapgene")
    if not rec.name or rec.name == "<unknown name>":
        rec.name = dna_path.stem            # SnapGene doesn't always store a name
    return rec


def read_trace(ab1_path: Path):
    """Read an .ab1 Sanger trace; carries per-base Phred quality."""
    return SeqIO.read(str(ab1_path), "abi")


def quality_trim(record) -> tuple[int, int]:
    """Return (start, end) of the high-quality region of the read.

    Trims the messy low-Q ends typical of Sanger reads without touching the
    interior, using a sliding-window mean and rejecting any window with an N.
    """
    q = record.letter_annotations["phred_quality"]
    seq = str(record.seq).upper()
    n = len(q)
    if n == 0:
        return 0, 0

    def window_ok(i: int) -> bool:
        w = q[i : i + TRIM_WINDOW]
        if not w or "N" in seq[i : i + TRIM_WINDOW]:
            return False
        return sum(w) / len(w) >= QUAL_TRIM_THRESHOLD

    start = 0
    while start < n and not window_ok(start):
        start += 1
    end = n
    while end > start and not window_ok(end - TRIM_WINDOW):
        end -= 1
    return start, max(start, end)


# --------------------------------------------------------------------------- #
# Alignment
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Discrepancy:
    ref_pos: int      # 1-based position in the reference (start, for indels)
    kind: str         # "mismatch" | "insertion" | "deletion"
    ref_base: str     # may be multi-base after consecutive deletions are collapsed
    read_base: str    # may be multi-base after consecutive insertions are collapsed
    read_qual: int    # Phred of the read base (or min flanking quality for a deletion)
    col: int          # column index in the gapped alignment (for context display)

    @property
    def confident(self) -> bool:
        """A real discrepancy vs. just an unread / low-quality base call."""
        if "N" in self.read_base.upper():
            return False
        return self.read_qual >= MIN_CONF_Q

    @property
    def span_len(self) -> int:
        """Number of bases involved (deleted, inserted, or 1 for a mismatch)."""
        if self.kind == "deletion":
            return len(self.ref_base)
        if self.kind == "insertion":
            return len(self.read_base)
        return 1

    @property
    def ref_range(self) -> str:
        """Human 1-based reference span, e.g. '3,031' or '4,000-4,002'."""
        if self.kind == "deletion" and self.span_len > 1:
            return f"{self.ref_pos:,}-{self.ref_pos + self.span_len - 1:,}"
        return f"{self.ref_pos:,}"


@dataclass(frozen=True)
class AlignmentResult:
    read_name: str
    ref_name: str
    strand: str               # "+" or "-"
    score: float
    ref_start: int            # 1-based, inclusive
    ref_end: int              # 1-based, inclusive
    aligned_length: int
    matches: int
    discrepancies: tuple[Discrepancy, ...]
    read_quality: tuple[int, ...]   # Phred over the aligned read
    aligned_ref: str                # gapped reference row (for context display)
    aligned_read: str               # gapped read row

    @property
    def confident_discrepancies(self) -> tuple[Discrepancy, ...]:
        # Filter to confident calls first, THEN collapse adjacent indels, so a
        # low-quality base between two real indels keeps them from merging.
        return _collapse_indels(tuple(d for d in self.discrepancies if d.confident))

    @property
    def percent_identity(self) -> float:
        return 100.0 * self.matches / self.aligned_length if self.aligned_length else 0.0

    @property
    def n_mismatch(self) -> int:
        return sum(1 for d in self.confident_discrepancies if d.kind == "mismatch")

    @property
    def n_indel(self) -> int:
        return sum(1 for d in self.confident_discrepancies if d.kind != "mismatch")

    @property
    def n_unread(self) -> int:
        return len(self.discrepancies) - len(self.confident_discrepancies)

    @property
    def is_perfect(self) -> bool:
        return not self.confident_discrepancies


def _collapse_indels(discreps: tuple[Discrepancy, ...]) -> tuple[Discrepancy, ...]:
    """Merge runs of consecutive single-base indels into one multi-base event.

    A 3 bp deletion arrives as three adjacent single-base deletion columns;
    this joins them into one "3 bp deletion". Adjacency is judged by alignment
    column so only a truly contiguous gap merges. Mismatches are never merged.
    """
    out: list[Discrepancy] = []
    for d in discreps:
        if out and d.kind == out[-1].kind and d.kind in ("insertion", "deletion"):
            prev = out[-1]
            prev_cols = len(prev.ref_base) if d.kind == "deletion" else len(prev.read_base)
            if d.col == prev.col + prev_cols:          # exactly the next column
                if d.kind == "deletion":
                    out[-1] = replace(prev, ref_base=prev.ref_base + d.ref_base,
                                      read_qual=min(prev.read_qual, d.read_qual))
                else:
                    out[-1] = replace(prev, read_base=prev.read_base + d.read_base,
                                      read_qual=min(prev.read_qual, d.read_qual))
                continue
        out.append(d)
    return tuple(out)


def _make_aligner() -> PairwiseAligner:
    aln = PairwiseAligner()
    aln.mode = "local"
    aln.match_score = MATCH
    aln.mismatch_score = MISMATCH
    aln.open_gap_score = OPEN_GAP
    aln.extend_gap_score = EXTEND_GAP
    return aln


def _walk_alignment(alignment, ref_offset, query_qual, query_offset):
    """Walk the gapped alignment columns; classify each into match or discrepancy.

    Emits 1-based reference coordinates so positions match SnapGene / IGV.
    """
    target = alignment[0]   # reference row, with '-' for gaps
    query = alignment[1]    # read row, with '-' for gaps
    ref_pos = ref_offset
    q_pos = query_offset
    found: list[Discrepancy] = []
    matches = 0
    for col, (t, q) in enumerate(zip(target, query)):
        if t == "-":                              # extra base in read (insertion)
            qv = query_qual[q_pos] if q_pos < len(query_qual) else 0
            found.append(Discrepancy(ref_pos, "insertion", "-", q, qv, col))
            q_pos += 1
        elif q == "-":                            # base missing from read (deletion)
            ref_pos += 1
            prev_q = query_qual[q_pos - 1] if 0 <= q_pos - 1 < len(query_qual) else 0
            next_q = query_qual[q_pos] if q_pos < len(query_qual) else 0
            found.append(Discrepancy(ref_pos, "deletion", t, "-", min(prev_q, next_q), col))
        else:
            ref_pos += 1
            qv = query_qual[q_pos] if q_pos < len(query_qual) else 0
            if t.upper() == q.upper():
                matches += 1
            else:
                found.append(Discrepancy(ref_pos, "mismatch", t, q, qv, col))
            q_pos += 1
    return tuple(found), matches


def align_read(ref, read) -> AlignmentResult:
    """Quality-trim, align both strands, and return the best AlignmentResult."""
    aligner = _make_aligner()
    start, end = quality_trim(read)
    trimmed = read.seq[start:end]
    trimmed_q = tuple(read.letter_annotations["phred_quality"][start:end])

    candidates = []
    # NOTE: reads align to the linear ref.seq (a REGION match), not a doubled
    # reference. A read that happens to span the map's coordinate origin would
    # align only to the larger side of the junction (partial coverage). This is
    # by design - Sanger reads cover a region; for whole-circle handling (an
    # assembly that wraps the origin) use compare_plasmid.py instead.
    for strand, seq in (("+", trimmed), ("-", trimmed.reverse_complement())):
        if len(seq) > 0:
            candidates.append((strand, aligner.align(ref.seq, seq)[0]))
    if not candidates:
        raise ValueError("read has no usable high-quality sequence after trimming")

    strand, aln = max(candidates, key=lambda c: c[1].score)
    quality = trimmed_q if strand == "+" else trimmed_q[::-1]
    ref_offset = int(aln.coordinates[0][0])
    query_offset = int(aln.coordinates[1][0])
    discrepancies, matches = _walk_alignment(aln, ref_offset, quality, query_offset)

    return AlignmentResult(
        read_name=read.id,
        ref_name=ref.name,
        strand=strand,
        score=aln.score,
        ref_start=ref_offset + 1,
        ref_end=int(aln.coordinates[0][-1]),
        aligned_length=matches + len(discrepancies),
        matches=matches,
        discrepancies=discrepancies,
        read_quality=quality,
        aligned_ref=str(aln[0]),
        aligned_read=str(aln[1]),
    )


# --------------------------------------------------------------------------- #
# Sequence-context display
# --------------------------------------------------------------------------- #
def context_lines(res: AlignmentResult, d: Discrepancy) -> tuple[str, str, str]:
    """Return (ref, match-bar, read) text around a discrepancy, gaps aligned."""
    lo = max(0, d.col - CONTEXT_BP)
    hi = d.col + CONTEXT_BP + 1
    ref_chunk = res.aligned_ref[lo:hi]
    read_chunk = res.aligned_read[lo:hi]
    bar = "".join(
        "|" if a != "-" and b != "-" and a.upper() == b.upper() else " "
        for a, b in zip(ref_chunk, read_chunk)
    )
    return ref_chunk, bar, read_chunk


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def print_report(res: AlignmentResult) -> None:
    conf = res.confident_discrepancies
    verdict = "PERFECT (clean match)" if res.is_perfect else f"{len(conf)} real discrepancy(ies)"
    print(f"\n{'='*70}")
    print(f"  {res.read_name}")
    print(f"  vs {res.ref_name}  (strand {res.strand})")
    print(f"{'='*70}")
    print(f"  Identity : {res.percent_identity:.2f}%  ({res.matches}/{res.aligned_length} aligned bases, raw)")
    print(f"  Covers   : reference {res.ref_start:,}-{res.ref_end:,}")
    print(f"  Mismatch : {res.n_mismatch}   Indels: {res.n_indel}   -> {verdict}")
    print(f"  ({res.n_unread} low-quality / N position(s) excluded from the verdict; identity is raw,")
    print(f"   so a PERFECT call can still read <100% when the read has low-quality bases)")
    for d in conf:
        label = {"mismatch": f"ref {d.ref_base} -> read {d.read_base}  [Q{d.read_qual}]",
                 "insertion": f"{d.span_len} bp inserted in read: {d.read_base}  [Q{d.read_qual}]",
                 "deletion": f"{d.span_len} bp missing from read (ref {d.ref_base})  [Q{d.read_qual}]"}[d.kind]
        ref_c, bar, read_c = context_lines(res, d)
        print(f"  {'-'*50}")
        print(f"  ref {d.ref_range}  {d.kind}: {label}")
        print(f"      ref  5'-{ref_c}-3'")
        print(f"             {bar}")
        print(f"      read 5'-{read_c}-3'")


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
_COLORS = {"mismatch": "#d62728", "insertion": "#1f77b4", "deletion": "#9467bd"}


def plot_read_detail(res: AlignmentResult, out_path: Path) -> None:
    """Per-read figure: discrepancy map + Phred track + sequence-context text."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    conf = res.confident_discrepancies
    fig, (ax_map, ax_q, ax_txt) = plt.subplots(
        3, 1, figsize=(12, 6.5), height_ratios=[1, 1.2, 1.6],
        gridspec_kw={"hspace": 0.5},
    )

    # Panel 1 - covered span with discrepancy ticks.
    ax_map.add_patch(plt.Rectangle((res.ref_start, 0.3), res.ref_end - res.ref_start,
                                   0.4, color="#2ca02c", alpha=0.25))
    ax_map.hlines(0.5, res.ref_start, res.ref_end, color="#2ca02c", lw=2)
    for d in conf:
        ax_map.vlines(d.ref_pos, 0.1, 0.9, color=_COLORS[d.kind], lw=1.5)
        if d.kind == "mismatch":
            tag = f"{d.ref_base}>{d.read_base}"
        else:
            tag = f"{d.kind[:3]}{d.span_len}" if d.span_len > 1 else d.kind[:3]
        ax_map.annotate(tag, (d.ref_pos, 0.93), fontsize=7, ha="center", color=_COLORS[d.kind])
    ax_map.set_xlim(res.ref_start - 20, res.ref_end + 20)
    ax_map.set_ylim(0, 1.05)
    ax_map.set_yticks([])
    ax_map.set_xlabel("reference position (bp)")
    seen = {d.kind for d in conf}
    if seen:
        ax_map.legend(handles=[Patch(color=_COLORS[k], label=k) for k in sorted(seen)],
                      loc="upper right", fontsize=8, ncol=3)

    # Panel 2 - Phred quality across the aligned read.
    q = res.read_quality
    ax_q.fill_between(range(len(q)), q, color="#888", alpha=0.5)
    ax_q.axhline(QUAL_TRIM_THRESHOLD, color="#d62728", ls="--", lw=0.8, label=f"Q{QUAL_TRIM_THRESHOLD}")
    ax_q.set_xlim(0, max(len(q), 1))
    ax_q.set_ylim(0, 65)
    ax_q.set_xlabel("read position (trimmed, bp)")
    ax_q.set_ylabel("Phred Q")
    ax_q.legend(loc="lower right", fontsize=8)

    # Panel 3 - actual aligned sequence around each discrepancy (monospace).
    ax_txt.axis("off")
    if res.is_perfect:
        ax_txt.text(0.5, 0.5, "Perfect match - no discrepancies in high-quality sequence",
                    ha="center", va="center", fontsize=11, color="#2ca02c", weight="bold")
    else:
        lines = []
        for d in conf[:6]:
            ref_c, bar, read_c = context_lines(res, d)
            head = (f"ref {d.ref_range}  {d.kind}"
                    + (f"  {d.ref_base}>{d.read_base} Q{d.read_qual}" if d.kind == "mismatch"
                       else f"  {d.span_len} bp"))
            lines += [head, f"  ref  {ref_c}", f"       {bar}", f"  read {read_c}", ""]
        if len(conf) > 6:
            lines.append(f"... and {len(conf) - 6} more (see console / summary)")
        ax_txt.text(0.01, 0.98, "\n".join(lines), ha="left", va="top",
                    fontfamily="monospace", fontsize=9, transform=ax_txt.transAxes)

    verdict = "PERFECT MATCH" if res.is_perfect else f"{res.n_mismatch} mismatch, {res.n_indel} indel"
    fig.suptitle(
        f"{res.read_name}  vs  {res.ref_name}  ({res.strand})\n"
        f"{res.percent_identity:.2f}% identity over ref {res.ref_start:,}-{res.ref_end:,}"
        f"   |   {verdict}",
        fontsize=10, y=0.99,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_combined(results: list[AlignmentResult], ref, out_path: Path) -> None:
    """One figure: every read tiled as a lane over the shared reference span."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    lo = min(r.ref_start for r in results)
    hi = max(r.ref_end for r in results)
    margin = max(20, (hi - lo) // 50)
    fig, ax = plt.subplots(figsize=(12, 1.0 + 0.5 * len(results)))

    for i, res in enumerate(results):
        y = len(results) - i
        ax.hlines(y, res.ref_start, res.ref_end, color="#2ca02c", lw=3, alpha=0.6)
        for d in res.confident_discrepancies:
            ax.plot(d.ref_pos, y, marker="v", color=_COLORS[d.kind], ms=6)
        tag = "OK" if res.is_perfect else f"{res.n_mismatch}mm/{res.n_indel}indel"
        ax.text(lo - margin, y, f"{res.read_name[:34]}  ({res.strand}, {tag})",
                ha="right", va="center", fontsize=8)

    ax.set_xlim(lo - margin, hi + margin)
    ax.set_ylim(0.3, len(results) + 0.7)
    ax.set_yticks([])
    ax.set_xlabel("reference position (bp)")
    ax.legend(handles=[Patch(color=_COLORS[k], label=k) for k in _COLORS],
              loc="upper right", fontsize=8, ncol=3)
    ax.set_title(f"{len(results)} reads vs {ref.name} ({len(ref.seq):,} bp)", fontsize=11)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _gather_reads(args: list[str]) -> list[Path]:
    reads: list[Path] = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            reads.extend(sorted(p.glob("*.ab1")))
        elif p.suffix.lower() == ".ab1" and p.exists():
            reads.append(p)
    return reads


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

    reads = _gather_reads(argv[1:])
    if not reads:
        print("ERROR: no .ab1 reads found. Pass .ab1 files or a folder containing them.")
        return 1

    out_dir = ref_path.parent / "alignment_figures"
    out_dir.mkdir(exist_ok=True)
    print(f"Reference: {ref.name} ({len(ref.seq):,} bp, {ref.annotations.get('topology', '?')})")
    print(f"Aligning {len(reads)} read(s); outputs -> {out_dir}")

    results: list[AlignmentResult] = []
    summary: list[str] = []
    for ab1 in reads:
        try:
            read = read_trace(ab1)
        except Exception as e:           # one bad read must not stop the batch
            print(f"\n  SKIPPED {ab1.name}: {e}")
            summary.append(f"SKIPPED   {ab1.name}: {e}")
            continue
        if len(read.seq) > SANGER_MAX_BP:
            print(f"\n  WARNING: {ab1.name} is {len(read.seq):,} bp - far longer than a Sanger read.")
            print(f"           If this is a WHOLE-PLASMID result (e.g. Plasmidsaurus), use")
            print(f"           compare_plasmid.py instead; this script only checks the part of")
            print(f"           the plasmid one read can span and may report a misleading match.")
        res = align_read(ref, read)
        print_report(res)
        plot_read_detail(res, out_dir / f"{_safe(res.read_name)}__vs__{_safe(ref.name)}.png")
        results.append(res)
        tag = "PERFECT" if res.is_perfect else f"{res.n_mismatch} mismatch, {res.n_indel} indel"
        summary.append(f"{tag:<22} {res.percent_identity:6.2f}%  ref {res.ref_start:,}-{res.ref_end:,}  {res.read_name}")

    if len(results) > 1:
        combined = out_dir / f"COMBINED__{_safe(ref.name)}.png"
        plot_combined(results, ref, combined)
        print(f"\nCombined figure -> {combined}")

    summary_path = out_dir / "alignment_summary.txt"
    header = f"Alignment summary  |  reference: {ref.name} ({len(ref.seq):,} bp)\n" + "=" * 70
    summary_path.write_text(header + "\n" + "\n".join(summary) + "\n", encoding="utf-8")
    print(f"Summary -> {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
