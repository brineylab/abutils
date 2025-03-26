#!usr/env/python
# filename: test_sequences.py

import os
import tempfile

import pandas as pd
import polars as pl
import pytest
from Bio.SeqRecord import SeqRecord

from ..core.sequence import Sequence
from ..io import (
    read_airr,
    read_csv,
    read_fasta,
    read_parquet,
    to_fasta,
)


@pytest.fixture
def fasta_sequence():
    f = ">fasta_sequence\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    fasta = tempfile.NamedTemporaryFile(delete=False)
    fasta.write(f.encode("utf-8"))
    fasta.close()
    return fasta


@pytest.fixture
def airr_sequences():
    a = "sequence_id\tsequence\tsequence_aa\nairr_sequence\tATCGATCGATCG\tCASSLDRYFDYW\nairr_sequence2\tATCGATCGATCG\tCASSLDRYFDYW\n"
    airr = tempfile.NamedTemporaryFile(delete=False)
    airr.write(a.encode("utf-8"))
    airr.close()
    return airr.name


@pytest.fixture
def parquet_sequences():
    a = [
        {
            "sequence_id": "airr_sequence",
            "sequence": "ATCGATCGATCG",
            "sequence_aa": "CASSLDRYFDYW",
        },
        {
            "sequence_id": "airr_sequence2",
            "sequence": "ATCGATCGATCG",
            "sequence_aa": "CASSLDRYFDYW",
        },
    ]
    df = pl.DataFrame(a)
    parquet = tempfile.NamedTemporaryFile(delete=False)
    df.write_parquet(parquet.name)
    parquet.close()
    return parquet.name


@pytest.fixture
def csv_sequences():
    a = "sequence_id,sequence,sequence_aa\nairr_sequence,ATCGATCGATCG,CASSLDRYFDYW\nairr_sequence2,ATCGATCGATCG,CASSLDRYFDYW\n"
    csv = tempfile.NamedTemporaryFile(delete=False)
    csv.write(a.encode("utf-8"))
    csv.close()
    return csv.name


@pytest.fixture
def polars_sequences():
    a = [
        {
            "sequence_id": "airr_sequence",
            "sequence": "ATCGATCGATCG",
            "sequence_aa": "CASSLDRYFDYW",
        },
        {
            "sequence_id": "airr_sequence2",
            "sequence": "ATCGATCGATCG",
            "sequence_aa": "CASSLDRYFDYW",
        },
    ]
    df = pl.DataFrame(a)
    return df


@pytest.fixture
def pandas_sequences():
    a = [
        {
            "sequence_id": "airr_sequence",
            "sequence": "ATCGATCGATCG",
            "sequence_aa": "CASSLDRYFDYW",
        },
        {
            "sequence_id": "airr_sequence2",
            "sequence": "ATCGATCGATCG",
            "sequence_aa": "CASSLDRYFDYW",
        },
    ]
    df = pd.DataFrame(a)
    return df


@pytest.fixture
def sequences():
    seq1 = Sequence(
        {
            "sequence": "ATCG",
            "seq": "ATCG",
            "sequence_id": "seq1",
            "name": "seq1",
        }
    )
    seq2 = Sequence(
        {
            "sequence": "GCTA",
            "seq": "GCTA",
            "sequence_id": "seq2",
            "name": "seq2",
        }
    )
    seq3 = Sequence(
        {
            "sequence": "TTTT",
            "seq": "TTTT",
            "sequence_id": "seq3",
            "name": "seq3",
        }
    )
    return [seq1, seq2, seq3]


@pytest.fixture
def seq_records():
    seq_records = [
        SeqRecord("ATCG", id="seq1"),
        SeqRecord("GCTA", id="seq2"),
        SeqRecord("TTTT", id="seq3"),
    ]
    return seq_records


# ------------------
#     Sequence
# ------------------


def test_sequence_with_string():
    seq = Sequence("ATCG")
    assert seq.sequence == "ATCG"
    assert seq.id is not None


def test_sequence_with_string_and_id():
    seq = Sequence("ATCG", id="seq1")
    assert seq.sequence == "ATCG"
    assert seq.id == "seq1"


def test_sequence_with_list():
    seq = Sequence(["seq1", "ATCG"])
    assert seq.sequence == "ATCG"
    assert seq.id == "seq1"


def test_sequence_with_dict():
    seq = Sequence({"sequence_id": "seq1", "sequence": "ATCG"})
    assert seq.sequence == "ATCG"
    assert seq.id == "seq1"


def test_sequence_with_dict_and_keys():
    seq = Sequence({"name": "seq1", "dna": "ATCG"}, id_key="name", seq_key="dna")
    assert seq.sequence == "ATCG"
    assert seq.id == "seq1"


def test_sequence_with_seqrecord():
    from Bio.SeqRecord import SeqRecord

    seq = SeqRecord("ATCG")
    seq.id = "seq1"
    ab_seq = Sequence(seq)
    assert ab_seq.sequence == "ATCG"
    assert ab_seq.id == "seq1"


def test_sequence_with_abutils_sequence():
    seq = Sequence("ATCG")
    ab_seq = Sequence(seq)
    assert ab_seq.sequence == "ATCG"
    assert ab_seq.id is not None


def test_sequence_contains():
    seq = Sequence("ATCG")
    assert "AT" in seq
    assert "CG" in seq
    assert "GC" not in seq


def test_sequence_getitem():
    seq = Sequence("ATCG")
    assert seq[0] == "A"
    assert seq[1:3] == "TC"
    assert seq[-1] == "G"


def test_sequence_fastq():
    seq = Sequence("ATCG", qual="abcd")
    assert seq.fastq == "@{}\n{}\n+\n{}".format(seq.id, seq.sequence, "abcd")


def test_sequence_reverse_complement():
    seq = Sequence("ATCG")
    assert seq.reverse_complement == "CGAT"


def test_sequence_annotations():
    seq = Sequence("ATCG")
    assert seq.annotations == {}
    seq.annotations = {"key": "value"}
    assert seq.annotations == {"key": "value"}


def test_sequence_strand():
    seq = Sequence("ATCG")
    assert seq.strand == "plus"
    seq.strand = "minus"
    assert seq.strand == "minus"


def test_sequence_translate():
    seq = Sequence("ATGTGCTAA")
    assert seq.translate() == "MC*"


def test_sequence_as_fasta():
    seq = Sequence("ATCG", id="seq1", annotations={"name": "my_seq", "seq": "TACG"})
    assert seq.as_fasta() == ">seq1\nATCG"
    assert seq.as_fasta(name_field="name") == ">my_seq\nATCG"
    assert seq.as_fasta(seq_field="seq") == ">seq1\nTACG"


def test_sequence_region():
    seq = Sequence("ATCG")
    assert seq.region() == ">{}\n{}".format(seq.id, seq.sequence)
    assert seq.region(start=1) == ">{}\n{}".format(seq.id, seq.sequence[1:])
    assert seq.region(end=2) == ">{}\n{}".format(seq.id, seq.sequence[:2])
    assert seq.region(start=1, end=2) == ">{}\n{}".format(seq.id, seq.sequence[1:2])


# ------------------
#     readers
# ------------------


def test_read_fasta(fasta_sequence):
    seqs = read_fasta(fasta_sequence.name)
    assert len(seqs) == 1
    assert seqs[0].id == "fasta_sequence"
    assert (
        seqs[0].sequence
        == "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    )


def test_read_airr(airr_sequences):
    seqs = read_airr(airr_sequences)
    assert len(seqs) == 2
    assert seqs[0].id == "airr_sequence"
    assert seqs[0].sequence == "ATCGATCGATCG"
    assert seqs[0]["sequence_aa"] == "CASSLDRYFDYW"
    assert seqs[1].id == "airr_sequence2"
    assert seqs[1].sequence == "ATCGATCGATCG"
    assert seqs[1]["sequence_aa"] == "CASSLDRYFDYW"


def test_read_parquet(parquet_sequences):
    seqs = read_parquet(parquet_sequences)
    assert len(seqs) == 2
    assert seqs[0].id == "airr_sequence"
    assert seqs[0].sequence == "ATCGATCGATCG"
    assert seqs[0]["sequence_aa"] == "CASSLDRYFDYW"
    assert seqs[1].id == "airr_sequence2"
    assert seqs[1].sequence == "ATCGATCGATCG"
    assert seqs[1]["sequence_aa"] == "CASSLDRYFDYW"


def test_read_csv(csv_sequences):
    seqs = read_csv(csv_sequences)
    assert len(seqs) == 2
    assert seqs[0].id == "airr_sequence"
    assert seqs[0].sequence == "ATCGATCGATCG"
    assert seqs[0]["sequence_aa"] == "CASSLDRYFDYW"
    assert seqs[1].id == "airr_sequence2"
    assert seqs[1].sequence == "ATCGATCGATCG"
    assert seqs[1]["sequence_aa"] == "CASSLDRYFDYW"


# def test_from_polars(polars_sequences):
#     seqs = from_polars(polars_sequences)
#     assert len(seqs) == 2
#     assert seqs[0].id == "airr_sequence"
#     assert seqs[0].sequence == "ATCGATCGATCG"
#     assert seqs[0]["sequence_aa"] == "CASSLDRYFDYW"
#     assert seqs[1].id == "airr_sequence2"
#     assert seqs[1].sequence == "ATCGATCGATCG"
#     assert seqs[1]["sequence_aa"] == "CASSLDRYFDYW"


# def test_from_pandas(pandas_sequences):
#     seqs = from_pandas(pandas_sequences)
#     assert len(seqs) == 2
#     assert seqs[0].id == "airr_sequence"
#     assert seqs[0].sequence == "ATCGATCGATCG"
#     assert seqs[0]["sequence_aa"] == "CASSLDRYFDYW"
#     assert seqs[1].id == "airr_sequence2"
#     assert seqs[1].sequence == "ATCGATCGATCG"
#     assert seqs[1]["sequence_aa"] == "CASSLDRYFDYW"


# ------------------
#    to fasta
# ------------------


def test_to_fasta_with_sequence_objects(sequences):
    fasta = to_fasta(sequences, as_string=True)
    assert isinstance(fasta, str)
    assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta


def test_to_fasta_with_fasta_string():
    fasta_string = ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT"
    fasta = to_fasta(fasta_string, as_string=True)
    assert isinstance(fasta, str)
    assert fasta == fasta_string


def test_to_fasta_with_fasta_file(sequences):
    fasta_string = "\n".join([s.fasta for s in sequences])
    fasta_file = tempfile.NamedTemporaryFile(delete=False)
    fasta_file.write(fasta_string.encode("utf-8"))
    fasta_file.close()
    fasta = to_fasta(fasta_file.name, as_string=True)
    assert isinstance(fasta, str)
    assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta


def test_to_fasta_with_seqrecord_objects(seq_records):
    fasta = to_fasta(seq_records, as_string=True)
    assert isinstance(fasta, str)
    assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta


def test_to_fasta_with_list_of_lists():
    seqs = [
        ["seq1", "ATCG"],
        ["seq2", "GCTA"],
        ["seq3", "TTTT"],
    ]
    fasta = to_fasta(seqs, as_string=True)
    assert isinstance(fasta, str)
    assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta


def test_to_fasta_with_id_key(sequences):
    fasta = to_fasta(sequences, id_key="name", as_string=True)
    assert isinstance(fasta, str)
    assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta


def test_to_fasta_with_sequence_key(sequences):
    fasta = to_fasta(sequences, sequence_key="seq", as_string=True)
    assert isinstance(fasta, str)
    assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta


def test_to_fasta_with_fasta_file_output(sequences):
    with tempfile.NamedTemporaryFile(delete=False) as fasta_file:
        fasta = to_fasta(sequences, fasta_file=fasta_file.name)
        assert isinstance(fasta, str)
        assert fasta == fasta_file.name
        with open(fasta_file.name, "r") as f:
            fasta_string = f.read()
        assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta_string


def test_to_fasta_with_fasta_string_output(sequences):
    fasta = to_fasta(sequences, as_string=True)
    assert isinstance(fasta, str)
    assert ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTTT" in fasta


# -------------------------------------------
#          codon optimization
# -------------------------------------------


def test_codon_optimize_returns_string():
    seq = Sequence("RINEYLA")
    optimized_seq = seq.codon_optimize(as_string=True)
    assert isinstance(optimized_seq, str)
    assert len(optimized_seq) == len(seq) * 3


def test_codon_optimize_returns_sequence_object():
    seq = Sequence("RINEYLA")
    optimized_seq = seq.codon_optimize(as_string=False)
    assert isinstance(optimized_seq, Sequence)
    assert len(optimized_seq.sequence) == len(seq) * 3
    assert optimized_seq.translate() == seq.sequence


def test_codon_optimize_with_frame():
    seq = Sequence("ATGCATGCATGCATGC")
    with pytest.raises(ValueError):
        seq.codon_optimize(frame=4)
