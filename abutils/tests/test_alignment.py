import pytest

from Bio.Align import MultipleSeqAlignment, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from abutils import Sequence
from abutils.tools.alignment import (
    consensus,
    mafft,
    muscle,
    PairwiseAlignment,
    LocalAlignment,
    GlobalAlignment,
    SemiGlobalAlignment,
    local_alignment,
    global_alignment,
    semiglobal_alignment,
    CIGAR,
    CIGARElement,
)


@pytest.fixture
def sequences():
    seq1 = SeqRecord(
        Seq("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"), id="seq1"
    )
    seq2 = SeqRecord(
        Seq("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"), id="seq2"
    )
    seq3 = SeqRecord(
        Seq("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"), id="seq3"
    )
    return [seq1, seq2, seq3]


@pytest.fixture
def aln_sequences():
    return MultipleSeqAlignment(
        [
            SeqRecord(
                Seq(
                    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
                    generic_protein,
                ),
                id="seq1",
            ),
            SeqRecord(
                Seq(
                    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
                    generic_protein,
                ),
                id="seq2",
            ),
            SeqRecord(
                Seq(
                    "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
                    generic_protein,
                ),
                id="seq3",
            ),
        ]
    )


# ---------------------------
#           MAFFT
# ---------------------------


def test_mafft_alignment(sequences):
    alignment = mafft(sequences)
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id == "seq1"
    assert alignment[1].id == "seq2"
    assert alignment[2].id == "seq3"


def test_mafft_alignment_file(sequences):
    alignment_file = "test_alignment.fasta"
    alignment = mafft(sequences, alignment_file=alignment_file, as_file=True)
    assert alignment == alignment_file


def test_mafft_alignment_string(sequences):
    alignment_string = mafft(sequences, as_string=True)
    assert isinstance(alignment_string, str)
    assert "seq1" in alignment_string
    assert "seq2" in alignment_string
    assert "seq3" in alignment_string


def test_mafft_alignment_phylip(sequences):
    alignment = mafft(sequences, fmt="phylip")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_mafft_alignment_clustal(sequences):
    alignment = mafft(sequences, fmt="clustal")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_mafft_alignment_threads(sequences):
    alignment = mafft(sequences, threads=2)
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_mafft_alignment_reorder(sequences):
    alignment = mafft(sequences, reorder=False)
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id == "seq1"
    assert alignment[1].id == "seq2"
    assert alignment[2].id == "seq3"


def test_mafft_alignment_id_key():
    seq1 = {
        "name": "seq1",
        "sequence": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq2 = {
        "name": "seq2",
        "sequence": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq3 = {
        "name": "seq3",
        "sequence": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    sequences = [seq1, seq2, seq3]
    alignment = mafft(sequences, id_key="name")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_mafft_alignment_seq_key():
    seq1 = {
        "id": "seq1",
        "seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq2 = {
        "id": "seq2",
        "seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq3 = {
        "id": "seq3",
        "seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    sequences = [seq1, seq2, seq3]
    alignment = mafft(sequences, seq_key="seq")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


# ---------------------------
#          MUSCLE
# ---------------------------


def test_muscle_alignment(sequences):
    alignment = muscle(sequences)
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id == "seq1"
    assert alignment[1].id == "seq2"
    assert alignment[2].id == "seq3"


def test_muscle_alignment_file(sequences):
    alignment_file = "test_alignment.fasta"
    alignment = muscle(sequences, alignment_file=alignment_file, as_file=True)
    assert alignment == alignment_file


def test_muscle_alignment_string(sequences):
    alignment_string = muscle(sequences, as_string=True)
    assert isinstance(alignment_string, str)
    assert "seq1" in alignment_string
    assert "seq2" in alignment_string
    assert "seq3" in alignment_string


def test_muscle_alignment_phylip(sequences):
    alignment = muscle(sequences, fmt="phylip")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_muscle_alignment_clustal(sequences):
    alignment = muscle(sequences, fmt="clustal")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_muscle_alignment_threads(sequences):
    alignment = muscle(sequences, threads=2)
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_muscle_alignment_reorder(sequences):
    alignment = muscle(sequences, reorder=False)
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id == "seq1"
    assert alignment[1].id == "seq2"
    assert alignment[2].id == "seq3"


def test_muscle_alignment_id_key():
    seq1 = {
        "name": "seq1",
        "sequence": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq2 = {
        "name": "seq2",
        "sequence": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq3 = {
        "name": "seq3",
        "sequence": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    sequences = [seq1, seq2, seq3]
    alignment = muscle(sequences, id_key="name")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


def test_muscle_alignment_seq_key():
    seq1 = {
        "id": "seq1",
        "seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq2 = {
        "id": "seq2",
        "seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    seq3 = {
        "id": "seq3",
        "seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
    }
    sequences = [seq1, seq2, seq3]
    alignment = muscle(sequences, seq_key="seq")
    assert isinstance(alignment, MultipleSeqAlignment)
    assert len(alignment) == 3
    assert alignment[0].id in ["seq1", "seq2", "seq3"]
    assert alignment[1].id in ["seq1", "seq2", "seq3"]
    assert alignment[2].id in ["seq1", "seq2", "seq3"]


# ---------------------------
#         CONSENSUS
# ---------------------------


def test_consensus(aln_sequences):
    name, consensus_string = consensus(aln_sequences)
    assert isinstance(name, str)
    assert isinstance(consensus_string, str)
    assert len(consensus_string) == 60
    assert (
        consensus_string
        == "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    )


# ---------------------------
#     PAIRWISE ALIGNMENTS
# ---------------------------


def test_pairwise_alignment_class():
    query = Sequence("ATCGATCGATCG")
    target = Sequence("ATCGATCGATCG")
    alignment = PairwiseAlignment(query, target, "test_target")
    assert alignment.query == query
    assert alignment.target == target
    assert alignment.target_id == "test_target"
    assert alignment.score is None


# ---------------------------
#    LOCAL ALIGNMENTS
# ---------------------------


def test_local_alignment_single_target():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = local_alignment(query, target=target)
    assert isinstance(alignment, LocalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_local_alignment_multiple_targets():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target1 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target2 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    targets = [target1, target2]
    alignments = local_alignment(query, targets=targets)
    assert isinstance(alignments, list)
    assert len(alignments) == 2
    for alignment in alignments:
        assert isinstance(alignment, LocalAlignment)
        assert isinstance(alignment.query, Sequence)
        assert isinstance(alignment.target, Sequence)
        assert len(alignment.aligned_query) == len(alignment.aligned_target)
        assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_local_alignment_seqrecord():
    query = SeqRecord("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = local_alignment(query, target=target)
    assert isinstance(alignment, LocalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_local_alignment_sequences():
    query = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    target = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    alignment = local_alignment(query, target)
    assert isinstance(alignment, LocalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_local_alignment_options():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = local_alignment(query, target=target, gap_open=-10, gap_extend=-0.5)
    assert isinstance(alignment, LocalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


# ---------------------------
#    GLOBAL ALIGNMENTS
# ---------------------------


def test_global_alignment_single_target():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = global_alignment(query, target=target)
    assert isinstance(alignment, GlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_global_alignment_multiple_targets():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target1 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target2 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    targets = [target1, target2]
    alignments = global_alignment(query, targets=targets)
    assert isinstance(alignments, list)
    assert len(alignments) == 2
    for alignment in alignments:
        assert isinstance(alignment, GlobalAlignment)
        assert isinstance(alignment.query, Sequence)
        assert isinstance(alignment.target, Sequence)
        assert len(alignment.aligned_query) == len(alignment.aligned_target)
        assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_global_alignment_seqrecord():
    query = SeqRecord("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = global_alignment(query, target=target)
    assert isinstance(alignment, GlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_global_alignment_sequences():
    query = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    target = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    alignment = global_alignment(query, target)
    assert isinstance(alignment, GlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_global_alignment_options():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = local_alignment(query, target=target, gap_open=-10, gap_extend=-0.5)
    assert isinstance(alignment, GlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


# ---------------------------
#   SEMIGLOBAL ALIGNMENTS
# ---------------------------


def test_semiglobal_alignment_single_target():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = semiglobal_alignment(query, target=target)
    assert isinstance(alignment, SemiGlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_semiglobal_alignment_multiple_targets():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target1 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target2 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    targets = [target1, target2]
    alignments = semiglobal_alignment(query, targets=targets)
    assert isinstance(alignments, list)
    assert len(alignments) == 2
    for alignment in alignments:
        assert isinstance(alignment, SemiGlobalAlignment)
        assert isinstance(alignment.query, Sequence)
        assert isinstance(alignment.target, Sequence)
        assert len(alignment.aligned_query) == len(alignment.aligned_target)
        assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_semiglobal_alignment_seqrecord():
    query = SeqRecord("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = semiglobal_alignment(query, target=target)
    assert isinstance(alignment, SemiGlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_semiglobal_alignment_sequences():
    query = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    target = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP")
    alignment = semiglobal_alignment(query, target)
    assert isinstance(alignment, SemiGlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


def test_semiglobal_alignment_options():
    query = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    target = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    alignment = local_alignment(query, target=target, gap_open=-10, gap_extend=-0.5)
    assert isinstance(alignment, SemiGlobalAlignment)
    assert isinstance(alignment.query, Sequence)
    assert isinstance(alignment.target, Sequence)
    assert len(alignment.aligned_query) == len(alignment.aligned_target)
    assert len(alignment.aligned_query) == len(alignment.alignment_midline)


# ---------------------------
#          CIGAR
# ---------------------------


def test_cigar_element():
    ce = CIGARElement(99, "=")
    assert ce.length == 99
    assert ce.element == "="


def test_cigar():
    cigar_string = "99=2X99="
    cigar = CIGAR(cigar_string)
    assert cigar.cigar_string == cigar_string
    cigar_list = cigar.cigar_list
    assert isinstance(cigar_list, list)
    assert len(cigar_list) == 4
    assert isinstance(cigar_list[0], CIGARElement)
    assert cigar_list[0].length == 99
    assert cigar_list[0].element == "="
    assert isinstance(cigar_list[1], CIGARElement)
    assert cigar_list[1].length == 2
    assert cigar_list[1].element == "X"
    assert isinstance(cigar_list[2], CIGARElement)
    assert cigar_list[2].length == 99
    assert cigar_list[2].element == "="
    assert isinstance(cigar_list[3], CIGARElement)
    assert cigar_list[3].length == 0
    assert cigar_list[3].element == ""
    assert len(cigar) == 200
    cigar_iter = list(cigar)
    assert isinstance(cigar_iter, list)
    assert len(cigar_iter) == 4
    assert isinstance(cigar_iter[0], CIGARElement)
    assert cigar_iter[0].length == 99
    assert cigar_iter[0].element == "="
    assert isinstance(cigar_iter[1], CIGARElement)
    assert cigar_iter[1].length == 2
    assert cigar_iter[1].element == "X"
    assert isinstance(cigar_iter[2], CIGARElement)
    assert cigar_iter[2].length == 99
    assert cigar_iter[2].element == "="
    assert isinstance(cigar_iter[3], CIGARElement)
    assert cigar_iter[3].length == 0
    assert cigar_iter[3].element == ""
    ce = cigar[1]
    assert isinstance(ce, CIGARElement)
    assert ce.length == 2
    assert ce.element == "X"
