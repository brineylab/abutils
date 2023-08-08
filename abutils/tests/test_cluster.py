import os
import pytest
import tempfile
import unittest

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
    seq1 = Sequence(seq=s1, id="seq1")
    seq2 = Sequence(seq=s1, id="seq2")
    seq3 = Sequence(seq=s1, id="seq3")
    seq4 = Sequence(seq=s2, id="seq4")
    seq5 = Sequence(seq=s2, id="seq5")
    return [seq1, seq2, seq3, seq4, seq5]


@pytest.fixture
def nt_sequences(nt_sequence_strings):
    # a total of 5 sequences, of which 3 are identical (s1) and 2 are identical (s2)
    s1, s2 = nt_sequence_strings
    seq1 = Sequence(seq=s1, id="seq1")
    seq2 = Sequence(seq=s1, id="seq2")
    seq3 = Sequence(seq=s1, id="seq3")
    seq4 = Sequence(seq=s2, id="seq4")
    seq5 = Sequence(seq=s2, id="seq5")
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
        SeqRecord(s1, "seq1"),
        SeqRecord(s1, "seq2"),
        SeqRecord(s1, "seq3"),
        SeqRecord(s2, "seq4"),
        SeqRecord(s2, "seq5"),
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
