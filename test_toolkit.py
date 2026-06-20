"""Regression tests for the plasmid-clone-validator toolkit.

Run from this folder:  python -m pytest test_toolkit.py -q

All test data is a synthetic, fabricated demo plasmid (example_test_data/), so
the suite is fully self-contained. It covers the documented known-answer cases
plus the load-bearing behaviors: completeness-gated verdicts, soft-clip
detection, circular (rotation) invariance, and minus-strand residue numbering.
"""
import re
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
DATA = HERE / "example_test_data"

import compare_plasmid as cp          # noqa: E402
import cloning_report as gr           # noqa: E402

REF = DATA / "demo_plasmid.gb"
CLEAN = DATA / "demo_clean.fasta"
MUTANT = DATA / "demo_mutant.fasta"
CAMPAIGN = DATA / "demo_campaign.txt"
CLONES = DATA / "demo_clones"

_COMPL = {"A": "C", "C": "A", "G": "T", "T": "G"}     # any non-identity base swap


def _query(seq_str, name):
    return SeqRecord(Seq(seq_str), id=name)


# --- documented known-answer cases ----------------------------------------- #
def test_clean_is_identical_and_fully_examined():
    res = cp.compare(cp.read_reference(REF), cp.read_query(CLEAN))
    assert res.is_identical and res.is_fully_examined
    assert res.clipped_query_bp == 0
    assert res.verdict == "IDENTICAL to map"


def test_mutant_known_answer():
    res = cp.compare(cp.read_reference(REF), cp.read_query(MUTANT))
    assert res.n_mismatch == 1 and res.n_indel == 1
    events = {(d.kind, d.ref_pos, d.span_len) for d in res.differences}
    assert ("mismatch", 1301, 1) in events           # planted SNP
    assert ("deletion", 1001, 3) in events            # planted 3 bp deletion


def test_gce_convergent_mutation():
    ref = cp.read_reference(REF)
    camp = gr.parse_config(CAMPAIGN)
    cds = gr.find_cds(ref, camp.gene)
    clones = [gr.analyze_clone(ref, cp.read_query(p), cds, camp)
              for p in sorted(CLONES.glob("*.fasta"))]
    conv = dict(gr.convergent_mutations(clones))
    assert conv.get("L25W", 0) == 3                   # planted convergent mutation


# --- completeness-gated verdict (the "earn the green light" guarantee) ------ #
def test_clipped_result_never_returns_clean_verdict():
    ref = cp.read_reference(REF)
    s = list(str(ref.seq).upper())
    s[-1] = _COMPL[s[-1]]                              # SNP at the very last base
    res = cp.compare(ref, _query("".join(s), "terminal"))
    assert res.clipped_query_bp >= 1
    assert res.is_identical and not res.is_fully_examined
    assert res.verdict != "IDENTICAL to map"
    assert "NOT a verified all-clear" in res.verdict


def test_interior_difference_not_clipped():
    ref = cp.read_reference(REF)
    s = list(str(ref.seq).upper())
    s[1500] = _COMPL[s[1500]]
    res = cp.compare(ref, _query("".join(s), "interior"))
    assert res.clipped_query_bp == 0 and res.n_mismatch == 1
    assert any(d.ref_pos == 1501 and d.kind == "mismatch" for d in res.differences)


def test_alignment_coordinates_full_circle():
    res = cp.compare(cp.read_reference(REF), cp.read_query(CLEAN))
    assert res.examined_ref_bp == res.ref_len and res.clipped_query_bp == 0


# --- circular handling ------------------------------------------------------ #
def test_rotation_invariance():
    ref = cp.read_reference(REF)
    s = list(str(ref.seq).upper())
    s[1500] = _COMPL[s[1500]]
    snp = "".join(s)
    unrot = cp.compare(ref, _query(snp, "unrot"))
    rot = cp.compare(ref, _query(snp[900:] + snp[:900], "rot"))   # cut the circle elsewhere
    a = sorted((d.kind, d.ref_pos) for d in unrot.differences)
    b = sorted((d.kind, d.ref_pos) for d in rot.differences)
    assert a == b == [("mismatch", 1501)]


# --- minus-strand CDS residue numbering ------------------------------------ #
def test_minus_strand_codon_residue():
    ref = cp.read_reference(REF)
    drug = next(f for f in ref.features
                if f.type == "CDS" and "drugr" in cp._flabel(f).lower())
    assert drug.location.strand == -1
    end = int(drug.location.end)
    pos0 = end - 1 - 3                                 # codon 2, first coding base
    fwd = str(ref.seq)[pos0].upper()
    con = cp._codon_consequence(ref, drug, pos0, _COMPL[fwd])
    assert con is not None and con.is_coding_change
    residue = int(re.search(r"p\.[A-Za-z*]+?(\d+)", con.text).group(1))
    assert residue == 2
    assert "DrugR" in con.text
