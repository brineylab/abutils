import os
import tempfile
import unittest

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..core.sequence import Sequence
from ..tools.cluster import (
    Cluster,
    Clusters,
    cluster,
    cluster_cdhit,
    cluster_mmseqs,
    cluster_vsearch,
    create_mmseqs_db,
    linclust,
)


@pytest.fixture
def aa_sequence_strings():
    s1 = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    s2 = "SADYMELTQKLTQEQLEDLRKCKENENTLVLVIRKPAKELQALSELKGDGQVNYDDDSSQ"
    return s1, s2


@pytest.fixture
def nt_sequence_strings():
    s1 = "ATCGGCTAATGCCGTAATGCCGTAATCGGCTAATGCCGTAATGCCGTAATCGGCTAATGCCGTAATGCCGTA"
    s2 = "CATGATTACTAGCGCTCATGATTACTAGCGCTCATGATTACTAGCGCTCATGATTACTAGCGCTCATGATTA"
    return s1, s2


@pytest.fixture
def aa_sequences(aa_sequence_strings):
    # a total of 5 sequences, of which 3 are identical (s1) and 2 are identical (s2)
    s1, s2 = aa_sequence_strings
    seq1 = Sequence(sequence=s1, id="seq1")
    seq2 = Sequence(sequence=s1, id="seq2")
    seq3 = Sequence(sequence=s1, id="seq3")
    seq4 = Sequence(sequence=s2, id="seq4")
    seq5 = Sequence(sequence=s2, id="seq5")
    return [seq1, seq2, seq3, seq4, seq5]


@pytest.fixture
def nt_sequences(nt_sequence_strings):
    # a total of 5 sequences, of which 3 are identical (s1) and 2 are identical (s2)
    s1, s2 = nt_sequence_strings
    seq1 = Sequence(sequence=s1, id="seq1")
    seq2 = Sequence(sequence=s1, id="seq2")
    seq3 = Sequence(sequence=s1, id="seq3")
    seq4 = Sequence(sequence=s2, id="seq4")
    seq5 = Sequence(sequence=s2, id="seq5")
    return [seq1, seq2, seq3, seq4, seq5]


@pytest.fixture
def aa_fasta_file(aa_sequence_strings):
    s1, s2 = aa_sequence_strings
    f = f">seq1\n{s1}\n"
    f += f">seq2\n{s1}\n"
    f += f">seq3\n{s1}\n"
    f += f">seq4\n{s2}\n"
    f += f">seq5\n{s2}\n"
    fasta = tempfile.NamedTemporaryFile(delete=False)
    fasta.write(f.encode("utf-8"))
    fasta.close()
    return fasta


@pytest.fixture
def nt_fasta_file(nt_sequence_strings):
    s1, s2 = nt_sequence_strings
    f = f">seq1\n{s1}\n"
    f += f">seq2\n{s1}\n"
    f += f">seq3\n{s1}\n"
    f += f">seq4\n{s2}\n"
    f += f">seq5\n{s2}\n"
    fasta = tempfile.NamedTemporaryFile(delete=False)
    fasta.write(f.encode("utf-8"))
    fasta.close()
    return fasta


def test_cluster_class(aa_sequences):
    seq1, seq2, seq3 = aa_sequences[:3]
    cluster = Cluster("cluster0", aa_sequences[:3])
    assert cluster.name == "cluster0"
    assert len(cluster.sequences) == 3
    assert seq1 in cluster.sequences
    assert seq2 in cluster.sequences
    assert seq3 in cluster.sequences
    assert cluster.size == 3
    assert cluster.seq_ids == ["seq1", "seq2", "seq3"]
    # consensus = cluster.consensus
    # assert consensus is not None
    # assert (
    #     str(consensus) == "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    # )


def test_clusters_class(aa_sequences):
    cluster1 = Cluster("cluster0", aa_sequences[:3])
    cluster2 = Cluster("cluster1", aa_sequences[3:])
    clusters = Clusters([cluster1, cluster2])
    assert len(clusters) == 2
    assert clusters.largest_cluster == cluster1
    centroids = clusters.centroids
    assert len(centroids) == 2
    assert cluster1.centroid in centroids
    assert cluster2.centroid in centroids
    assert clusters.clusters == [cluster1, cluster2]
    for i, cluster in enumerate(clusters):
        assert cluster == clusters[i]
    assert clusters[0] == cluster1
    assert clusters[1] == cluster2


def test_cluster_mmseqs(aa_fasta_file):
    cluster_info = cluster_mmseqs(
        aa_fasta_file.name,
        threshold=0.99,
        temp_dir=tempfile.mkdtemp(),
        as_dict=True,
        debug=True,
    )
    assert len(cluster_info) == 2
    name = list(cluster_info.keys())[0]
    assert len(cluster_info[name]["seq_ids"]) == 3
    os.unlink(aa_fasta_file.name)


def test_cluster_cdhit(aa_fasta_file):
    cluster_info = cluster_cdhit(
        aa_fasta_file.name,
        threshold=0.99,
        temp_dir=tempfile.mkdtemp(),
        as_dict=True,
        debug=True,
    )
    assert len(cluster_info) == 2
    name = list(cluster_info.keys())[0]
    assert len(cluster_info[name]["seq_ids"]) == 3
    os.unlink(aa_fasta_file.name)


def test_cluster_vsearch(nt_fasta_file):
    cluster_info = cluster_vsearch(
        nt_fasta_file.name,
        threshold=0.99,
        temp_dir=tempfile.mkdtemp(),
        as_dict=True,
        debug=True,
    )
    assert len(cluster_info) == 2
    name = list(cluster_info.keys())[0]
    assert len(cluster_info[name]["seq_ids"]) == 3
    os.unlink(nt_fasta_file.name)


def test_cluster_with_list_of_sequences(nt_sequence_strings):
    s1, s2 = nt_sequence_strings
    seqs = [s1, s1, s1, s2, s2]
    clusters = cluster(seqs, threshold=0.9)
    assert isinstance(clusters, Clusters)
    assert len(clusters) == 2
    assert clusters.count == 2
    assert clusters.clusters[0].size == 3
    assert clusters.clusters[1].size == 2


def test_cluster_with_fasta_file(nt_fasta_file):
    clusters = cluster(nt_fasta_file.name, threshold=0.9)
    assert isinstance(clusters, Clusters)
    assert len(clusters) == 2
    assert clusters.count == 2
    assert clusters.clusters[0].size == 3
    assert clusters.clusters[1].size == 2
    os.unlink(nt_fasta_file.name)


def test_cluster_with_seqrecord_objects(nt_sequence_strings):
    s1, s2 = nt_sequence_strings
    seqs = [
        SeqRecord(Seq(s1), id="seq1"),
        SeqRecord(Seq(s1), id="seq2"),
        SeqRecord(Seq(s1), id="seq3"),
        SeqRecord(Seq(s2), id="seq4"),
        SeqRecord(Seq(s2), id="seq5"),
    ]
    clusters = cluster(seqs, threshold=0.9)
    assert isinstance(clusters, Clusters)
    assert len(clusters) == 2
    assert clusters.count == 2
    assert clusters.clusters[0].size == 3
    assert clusters.clusters[1].size == 2


def test_cluster_with_sequence_objects(nt_sequences):
    clusters = cluster(nt_sequences, threshold=0.9)
    assert isinstance(clusters, Clusters)
    assert len(clusters) == 2
    assert clusters.count == 2
    assert clusters.clusters[0].size == 3
    assert clusters.clusters[1].size == 2


def test_cluster_with_list_of_lists(nt_sequence_strings):
    s1, s2 = nt_sequence_strings
    seqs = [["seq1", s1], ["seq2", s1], ["seq3", s1], ["seq4", s2], ["seq5", s2]]
    clusters = cluster(seqs, threshold=0.9)
    assert isinstance(clusters, Clusters)
    assert len(clusters) == 2
    assert clusters.count == 2
    assert clusters.clusters[0].size == 3
    assert clusters.clusters[1].size == 2


# -----------------
# MMseqs2 database and linclust tests
# -----------------


def test_create_mmseqs_db(aa_fasta_file):
    """Test database creation from FASTA file."""
    temp_dir = tempfile.mkdtemp()
    db_path = os.path.join(temp_dir, "test_db")
    result = create_mmseqs_db(aa_fasta_file.name, db_path, quiet=True)
    assert result == db_path
    assert os.path.exists(db_path)
    assert os.path.exists(f"{db_path}.index")
    assert os.path.exists(f"{db_path}.dbtype")
    os.unlink(aa_fasta_file.name)


def test_create_mmseqs_db_from_sequences(aa_sequences):
    """Test database creation from Sequence objects."""
    temp_dir = tempfile.mkdtemp()
    db_path = os.path.join(temp_dir, "test_db")
    result = create_mmseqs_db(aa_sequences, db_path, quiet=True)
    assert result == db_path
    assert os.path.exists(db_path)
    assert os.path.exists(f"{db_path}.index")


def test_linclust_basic(aa_fasta_file):
    """Test basic linclust functionality with TSV output."""
    temp_dir = tempfile.mkdtemp()
    output_tsv = os.path.join(temp_dir, "clusters.tsv")

    result = linclust(
        aa_fasta_file.name,
        output_tsv=output_tsv,
        threshold=0.9,
        quiet=True,
    )

    assert result["tsv_file"] == output_tsv
    assert result["fasta_file"] is None
    assert os.path.exists(output_tsv)
    assert result["cluster_count"] == 2
    os.unlink(aa_fasta_file.name)


def test_linclust_with_fasta_output(aa_fasta_file):
    """Test linclust with representative FASTA output."""
    temp_dir = tempfile.mkdtemp()
    output_tsv = os.path.join(temp_dir, "clusters.tsv")
    output_fasta = os.path.join(temp_dir, "reps.fasta")

    result = linclust(
        aa_fasta_file.name,
        output_tsv=output_tsv,
        output_fasta=output_fasta,
        threshold=0.9,
        quiet=True,
    )

    assert result["tsv_file"] == output_tsv
    assert result["fasta_file"] == output_fasta
    assert os.path.exists(output_tsv)
    assert os.path.exists(output_fasta)
    assert result["cluster_count"] == 2

    # Verify FASTA file has the expected number of sequences (2 representatives)
    with open(output_fasta) as f:
        fasta_content = f.read()
    rep_count = fasta_content.count(">")
    assert rep_count == 2
    os.unlink(aa_fasta_file.name)


def test_linclust_with_pregenerated_db(aa_sequences):
    """Test linclust with pre-generated database."""
    temp_dir = tempfile.mkdtemp()
    db_path = os.path.join(temp_dir, "test_db")
    output_tsv = os.path.join(temp_dir, "clusters.tsv")

    # Create database first
    create_mmseqs_db(aa_sequences, db_path, quiet=True)

    # Run linclust with the pre-generated database
    result = linclust(db_path, output_tsv=output_tsv, threshold=0.9, quiet=True)

    assert os.path.exists(output_tsv)
    assert result["cluster_count"] == 2


def test_linclust_tsv_format(aa_fasta_file):
    """Test that TSV output is in native mmseqs format (two columns)."""
    temp_dir = tempfile.mkdtemp()
    output_tsv = os.path.join(temp_dir, "clusters.tsv")

    linclust(aa_fasta_file.name, output_tsv=output_tsv, threshold=0.9, quiet=True)

    with open(output_tsv) as f:
        lines = f.readlines()

    # Verify format: two tab-separated columns
    assert len(lines) > 0
    for line in lines:
        if line.strip() and not line.startswith("#"):
            parts = line.strip().split("\t")
            assert len(parts) == 2, "TSV should have exactly 2 columns"
    os.unlink(aa_fasta_file.name)


def test_cluster_mmseqs_with_db_path(aa_sequences):
    """Test cluster_mmseqs with pre-generated database."""
    temp_dir = tempfile.mkdtemp()
    db_path = os.path.join(temp_dir, "test_db")

    # Create database first
    create_mmseqs_db(aa_sequences, db_path, quiet=True)

    # Run cluster_mmseqs with the pre-generated database
    cluster_info = cluster_mmseqs(
        db_path=db_path,
        threshold=0.99,
        temp_dir=temp_dir,
        as_dict=True,
        debug=False,
    )

    assert len(cluster_info) == 2


def test_cluster_mmseqs_validation_error():
    """Test that cluster_mmseqs raises ValueError when neither input is provided."""
    with pytest.raises(ValueError) as excinfo:
        cluster_mmseqs(threshold=0.99)
    assert "Either 'fasta_file' or 'db_path' must be provided" in str(excinfo.value)


def test_cluster_mmseqs_both_inputs_error(aa_fasta_file, aa_sequences):
    """Test that cluster_mmseqs raises ValueError when both inputs are provided."""
    temp_dir = tempfile.mkdtemp()
    db_path = os.path.join(temp_dir, "test_db")
    create_mmseqs_db(aa_sequences, db_path, quiet=True)

    with pytest.raises(ValueError) as excinfo:
        cluster_mmseqs(fasta_file=aa_fasta_file.name, db_path=db_path, threshold=0.99)
    assert "Only one of 'fasta_file' or 'db_path' should be provided" in str(
        excinfo.value
    )
    os.unlink(aa_fasta_file.name)
