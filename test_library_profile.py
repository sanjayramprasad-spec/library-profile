"""Regression tests for library_profile.py.

Run from this folder:  python -m pytest test_chimera.py -q

All data is synthetic and built in-memory: a tiny panel of source references that
share a backbone + two conserved junctions and differ inside three domains, plus
reads stitched from known domain mosaics (with planted Nanopore-style errors). The
suite checks the load-bearing behaviors: domain auto-detection, known-answer
genotype calls, strand invariance, and the honesty gate (partial / ambiguous /
unassigned reads never count toward composition).
"""
import random
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import library_profile as cpf          # noqa: E402

BASES = "ACGT"
_COMPL = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMPL)[::-1]


def _rand(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(BASES) for _ in range(n))


# A backbone, two conserved junctions, and 3 variable domains x 4 sources.
BACKBONE5 = _rand(120, 1)
BACKBONE3 = _rand(400, 2)
JUNCTION1 = _rand(30, 3)           # conserved N<->cat junction
JUNCTION2 = _rand(30, 4)           # conserved cat<->C junction
SOURCES = ["srcA", "srcB", "srcC", "srcD"]

# Each source's three domain alleles are independent random sequences, so every
# allele is ~75% divergent from the others -> rich source-specific k-mers.
DOMAINS = {"N": 240, "cat": 180, "C": 300}
ALLELES = {
    d: {s: _rand(L, 100 + 7 * di + 3 * si) for si, s in enumerate(SOURCES)}
    for di, (d, L) in enumerate(DOMAINS.items())
}


def _ref(source: str) -> SeqRecord:
    seq = (BACKBONE5 + ALLELES["N"][source] + JUNCTION1
           + ALLELES["cat"][source] + JUNCTION2 + ALLELES["C"][source] + BACKBONE3)
    return SeqRecord(Seq(seq), id=source, description="")


def _refs() -> list:
    return [_ref(s) for s in SOURCES]


def _chimera(n_src: str, cat_src: str, c_src: str, *, err: float = 0.0,
             seed: int = 0, full: bool = True) -> str:
    seq = (BACKBONE5 + ALLELES["N"][n_src] + JUNCTION1 + ALLELES["cat"][cat_src]
           + JUNCTION2 + ALLELES["C"][c_src] + BACKBONE3)
    if not full:                       # truncate to N + catalytic only
        seq = BACKBONE5 + ALLELES["N"][n_src] + JUNCTION1 + ALLELES["cat"][cat_src]
    if err:
        rng = random.Random(seed)
        chars = list(seq)
        for i in range(len(chars)):
            if rng.random() < err:
                chars[i] = rng.choice(BASES)
        seq = "".join(chars)
    return seq


def _panel():
    return cpf.build_panel(_refs(), k=15, anchor_min=20, names=("N", "cat", "C"))


# --- domain auto-detection -------------------------------------------------- #
def test_detects_three_domains_between_junctions():
    panel = _panel()
    assert [d.name for d in panel.domains] == ["N", "cat", "C"]
    # every source has private markers in every domain
    for d in panel.domains:
        assert all(panel.marker_pool[d.name][s] > 0 for s in SOURCES)


def test_names_count_must_match_detected_domains():
    import pytest
    with pytest.raises(ValueError):
        cpf.build_panel(_refs(), k=15, names=("only", "two"))


# --- known-answer genotype calls ------------------------------------------- #
def test_clean_chimera_known_answer():
    panel = _panel()
    seq = _chimera("srcA", "srcC", "srcB")
    call = cpf.call_read("r1", seq, panel)
    assert call.status == "complete"
    assert call.genotype == ("srcA", "srcC", "srcB")


def test_noisy_chimera_still_called():
    panel = _panel()
    seq = _chimera("srcD", "srcA", "srcC", err=0.05, seed=42)
    call = cpf.call_read("r2", seq, panel)
    assert call.status == "complete"
    assert call.genotype == ("srcD", "srcA", "srcC")


def test_strand_invariance():
    panel = _panel()
    seq = _chimera("srcB", "srcD", "srcA")
    fwd = cpf.call_read("f", seq, panel)
    rev = cpf.call_read("r", _rc(seq), panel)
    assert fwd.genotype == rev.genotype == ("srcB", "srcD", "srcA")
    assert fwd.strand == "+" and rev.strand == "-"


# --- honesty gate ----------------------------------------------------------- #
def test_truncated_read_is_partial_not_counted():
    panel = _panel()
    seq = _chimera("srcA", "srcB", "srcC", full=False)   # missing the C domain
    call = cpf.call_read("t", seq, panel)
    assert call.status == "partial"
    assert call.genotype is None
    assert call.calls["C"].source is None                # the absent domain
    assert call.calls["N"].source == "srcA"              # present domains still resolved


def test_unknown_source_domain_flagged_partial():
    panel = _panel()
    # A C-term from a source NOT in the panel -> no markers hit -> unassigned.
    novel_c = _rand(300, 9999)
    seq = (BACKBONE5 + ALLELES["N"]["srcA"] + JUNCTION1 + ALLELES["cat"]["srcB"]
           + JUNCTION2 + novel_c + BACKBONE3)
    call = cpf.call_read("u", seq, panel)
    assert call.calls["C"].source is None
    assert call.status in ("partial", "ambiguous")
    assert call.genotype is None


def test_composition_excludes_incomplete_reads():
    panel = _panel()
    recs = [
        SeqRecord(Seq(_chimera("srcA", "srcB", "srcC")), id="g1"),
        SeqRecord(Seq(_chimera("srcA", "srcB", "srcC")), id="g2"),
        SeqRecord(Seq(_chimera("srcD", "srcD", "srcD")), id="g3"),
        SeqRecord(Seq(_chimera("srcA", "srcB", "srcC", full=False)), id="partial"),
    ]
    res = cpf.profile_reads(panel, recs)
    assert len(res.complete) == 3
    comp = dict(res.composition())
    assert comp[("srcA", "srcB", "srcC")] == 2
    assert comp[("srcD", "srcD", "srcD")] == 1
    assert res.status_counts()["partial"] == 1


def test_domain_usage_counts():
    panel = _panel()
    recs = [SeqRecord(Seq(_chimera("srcA", "srcB", "srcC")), id=f"r{i}") for i in range(5)]
    res = cpf.profile_reads(panel, recs)
    usage = res.domain_usage()
    assert usage["N"]["srcA"] == 5
    assert usage["cat"]["srcB"] == 5
    assert usage["C"]["srcC"] == 5


def _blended_cat(a: str, b: str) -> str:
    """A read whose catalytic domain is half source a, half source b -> roughly
    equal private-marker hits for both, which must trip the margin (ambiguous) gate."""
    cat = ALLELES["cat"][a][: DOMAINS["cat"] // 2] + ALLELES["cat"][b][DOMAINS["cat"] // 2:]
    return (BACKBONE5 + ALLELES["N"]["srcA"] + JUNCTION1 + cat
            + JUNCTION2 + ALLELES["C"]["srcC"] + BACKBONE3)


def test_close_margin_is_ambiguous_and_excluded():
    panel = _panel()
    call = cpf.call_read("amb", _blended_cat("srcA", "srcB"), panel)
    assert call.calls["cat"].source is None          # neither source cleared the margin
    assert call.status == "ambiguous"
    assert call.genotype is None
    res = cpf.profile_reads(panel, [SeqRecord(Seq(_blended_cat("srcA", "srcB")), id="x")])
    assert len(res.complete) == 0
    assert res.status_counts()["ambiguous"] == 1


def test_partial_outranks_ambiguous_in_status():
    panel = _panel()
    # catalytic blended (ambiguous) AND C-term from a source not in the panel (partial).
    cat = ALLELES["cat"]["srcA"][:90] + ALLELES["cat"]["srcB"][90:]
    novel_c = _rand(300, 7777)
    seq = (BACKBONE5 + ALLELES["N"]["srcA"] + JUNCTION1 + cat
           + JUNCTION2 + novel_c + BACKBONE3)
    call = cpf.call_read("mix", seq, panel)
    assert call.calls["cat"].source is None and call.calls["C"].source is None
    assert call.status == "partial"                  # partial wins over ambiguous
    assert call.genotype is None


def test_names_parsing_accepts_comma_and_space_forms():
    # All the ways a shell can hand us three names must yield the same tuple.
    base = ["sources.fasta", "reads.fastq", "--names"]
    expected = ("N-term", "catalytic", "C-term")
    assert cpf._parse_args(base + ["N-term,catalytic,C-term"])["names"] == expected
    assert cpf._parse_args(base + ["N-term", "catalytic", "C-term"])["names"] == expected
    assert cpf._parse_args(base + ["N-term,", "catalytic,", "C-term"])["names"] == expected
    # --names must not swallow the positional read/ref paths.
    parsed = cpf._parse_args(["refs.fasta", "--names", "N", "mid", "C", "reads.fastq"])
    assert parsed["names"] == ("N", "mid", "C")
    assert parsed["positional"] == ["refs.fasta", "reads.fastq"]


def test_missing_option_value_errors_cleanly():
    import pytest
    with pytest.raises(SystemExit):
        cpf._parse_args(["refs.fasta", "reads.fastq", "--k"])


def test_thin_marker_sources_are_warned():
    panel = _panel()
    # Real pools are a few hundred markers; a normal floor flags nothing...
    assert cpf.thin_marker_warnings(panel, min_markers=10) == []
    # ...but an absurd floor flags every source in every domain (the silent-bias guard).
    warns = cpf.thin_marker_warnings(panel, min_markers=10_000)
    assert len(warns) == len(panel.sources) * len(panel.domains)
