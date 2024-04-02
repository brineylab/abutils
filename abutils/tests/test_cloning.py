import pytest

from abutils import Sequence
from abutils.tools.cloning import build_synthesis_constructs

# ------------------------------------
#      build_synthesis_constructs
# ------------------------------------


def test_build_synthesis_constructs_single_sequence():
    seq = Sequence("RINEYLA")
    seq["locus"] = "IGH"
    constructs = build_synthesis_constructs(seq)
    assert len(constructs) == 1
    opt_seq = str(constructs[0].sequence)
    assert opt_seq.startswith("catcctttttctagtagcaactgcaaccggtgtacac".upper())
    assert opt_seq.endswith("gcgtcgaccaagggcccatcggtcttcc".upper())


def test_build_synthesis_constructs_multiple_sequences():
    seqs = [Sequence("RINEYLA"), Sequence("ALYENIR")]
    for seq in seqs:
        seq["locus"] = "IGH"
    constructs = build_synthesis_constructs(seqs)
    assert len(constructs) == 2


def test_build_synthesis_constructs_with_overhangs():
    seq = Sequence("RINEYLA")
    seq["locus"] = "IGH"
    overhang_5 = {"IGH": "gggg"}
    overhang_3 = {"IGH": "cccc"}
    constructs = build_synthesis_constructs(
        seq, overhang_5=overhang_5, overhang_3=overhang_3
    )
    opt_seq = str(constructs[0].sequence)
    assert opt_seq.startswith("gggg".upper())
    assert opt_seq.endswith("cccc".upper())


def test_build_synthesis_constructs_no_overhangs():
    seq = Sequence("ATGCATGCTTATGAG")
    constructs = build_synthesis_constructs(seq, overhang_5={}, overhang_3={})
    assert len(constructs[0].sequence) == len("ATGCATGCTTATGAG")


def test_build_synthesis_constructs_group_by_chain():
    seqs = [Sequence("RINEYLA"), Sequence("ALYENIR")]
    for seq, locus in zip(seqs, ["IGH", "IGK"]):
        seq["locus"] = locus
    constructs = build_synthesis_constructs(seqs, group_by_chain=True)
    assert constructs[0].id[-3:] == "IGH"
    assert constructs[1].id[-3:] == "IGK"


def test_build_synthesis_constructs_sort_func():
    seqs = [Sequence("RINEYLA", id="B"), Sequence("ALYENIR", id="A")]
    for s in seqs:
        s["locus"] = "IGH"
    constructs = build_synthesis_constructs(
        seqs, sort_func=lambda s: s.id, add_locus_to_name=False
    )
    assert constructs[0].id == "A"
    assert constructs[1].id == "B"
