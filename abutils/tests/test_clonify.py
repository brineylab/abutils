import pytest

from abutils import Sequence
from abutils.tools.clonify import clonify, pairwise_distance


@pytest.fixture
def clonify_seqs():
    return [
        Sequence(
            id="seq1",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq2",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq3",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq4",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq5",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq6",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq7",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq8",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq9",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
        Sequence(
            id="seq10",
            seq="ATCGATCGATCG",
            v_gene="IGHV1-2",
            j_gene="IGHJ4",
            cdr3_aa="CASSLDRYFDYW",
        ),
    ]


@pytest.fixture
def pairwise_distance_seqs():
    seq1 = Sequence(
        "seq1", "ATCGATCGATCG", cdr3="CASSLDRYFDYW", mutations=["V1A", "V2G", "V3C"]
    )
    seq2 = Sequence(
        "seq2", "ATCGATCGATCG", cdr3="CASSLDRYFDYW", mutations=["V1A", "V2G", "V3C"]
    )
    seq3 = Sequence(
        "seq3", "ATCGATCGATCG", cdr3="CASSLDYFDYW", mutations=["V1A", "V2G", "V3T"]
    )
    return seq1, seq2, seq3


def test_clonify(clonify_seqs):
    results = clonify(clonify_seqs)
    assert isinstance(results, dict)
    assert len(results) == 1
    lineage = list(results.values())[0]
    assert isinstance(lineage, list)
    assert len(lineage) == len(clonify_seqs)
    for seq in clonify_seqs:
        assert seq in lineage


def test_pairwise_distance(pairwise_distance_seqs):
    seq1, seq2, seq3 = pairwise_distance_seqs
    dist1 = pairwise_distance(seq1, seq2)
    dist2 = pairwise_distance(seq1, seq3)
    assert isinstance(dist1, float)
    assert isinstance(dist2, float)
    assert dist1 <= 0.01
    assert dist2 > 0.0
