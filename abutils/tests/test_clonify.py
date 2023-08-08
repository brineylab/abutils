from collections import Counter
import pytest

from ..core.sequence import Sequence
from ..tools.clonify import clonify, pairwise_distance


@pytest.fixture
def clonify_seqs():
    return [
        Sequence(
            dict(
                sequence_id="seq1",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq2",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq3",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq4",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq5",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq6",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq7",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq8",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-2",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq9",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-69",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
        Sequence(
            dict(
                sequence_id="seq10",
                sequence="ATCGATCGATCG",
                v_gene="IGHV1-69",
                j_gene="IGHJ4",
                cdr3_aa="CASSLDRYFDYW",
                v_mutations=[],
                locus="IGH",
            )
        ),
    ]


@pytest.fixture
def pairwise_distance_seqs():
    seq1 = Sequence(
        dict(
            sequence_id="seq1",
            sequence="ATCGATCGATCG",
            cdr3="CASSLDRYFDYW",
            mutations=["V1A", "V2G", "V3C"],
        )
    )
    seq2 = Sequence(
        dict(
            sequence_id="seq2",
            sequence="ATCGATCGATCG",
            cdr3="CASSLDRYFDYW",
            mutations=["V1A", "V2G", "V3C"],
        )
    )
    seq3 = Sequence(
        dict(
            sequence_id="seq3",
            sequence="ATCGATCGATCG",
            cdr3="CASSLDYFDYW",
            mutations=["V1A", "V2G", "V3T"],
        )
    )
    return seq1, seq2, seq3


def test_clonify(clonify_seqs):
    results = clonify(clonify_seqs, return_assignment_dict=True)
    # sanity checks on the results
    assert isinstance(results, dict)
    for seq in clonify_seqs:
        assert seq.id in results
    # check for correct lineage sizes
    lineage_sizes = Counter(results.values())
    lineage_sizes = dict(
        sorted(lineage_sizes.items(), key=lambda x: x[1], reverse=True)
    )
    assert len(lineage_sizes) == 2
    lineage_names = list(lineage_sizes.keys())
    assert lineage_sizes[lineage_names[0]] == 8
    assert lineage_sizes[lineage_names[1]] == 2


def test_pairwise_distance(pairwise_distance_seqs):
    seq1, seq2, seq3 = pairwise_distance_seqs
    dist1 = pairwise_distance(seq1, seq2)
    dist2 = pairwise_distance(seq1, seq3)
    assert isinstance(dist1, float)
    assert isinstance(dist2, float)
    assert dist1 <= 0.01
    assert dist2 > 0.0
