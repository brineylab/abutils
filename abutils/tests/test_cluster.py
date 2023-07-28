import os
import pytest
import tempfile
import unittest

from Bio.Seq import Seq, SeqRecord

from abutils import Sequence
from abutils.tools.cluster import Cluster, Clusters
from abutils.tl import cluster, cluster_mmseqs, cluster_vsearch


@pytest.fixture
def sequences():
    seq1 = Sequence(
        seq="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
        id="seq1",
    )
    seq2 = Sequence(
        seq="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
        id="seq2",
    )
    seq3 = Sequence(
        seq="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
        id="seq3",
    )
    seq4 = Sequence(
        seq="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
        id="seq4",
    )
    seq5 = Sequence(
        seq="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP",
        id="seq5",
    )
    return [seq1, seq2, seq3, seq4, seq5]


@pytest.fixture
def fasta_file():
    f = ">seq1\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
    f += ">seq2\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
    f += ">seq3\nMEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
    fasta = tempfile.NamedTemporaryFile(delete=False)
    fasta.write(f.encode("utf-8"))
    fasta.close()
    return fasta


def test_cluster_class(sequences):
    seq1, seq2, seq3 = sequences[:3]
    cluster = Cluster("cluster0", sequences)
    assert cluster.name == "cluster0"
    assert len(cluster.sequences) == 3
    assert seq1 in cluster.sequences
    assert seq2 in cluster.sequences
    assert seq3 in cluster.sequences
    assert cluster.size == 3
    assert cluster.seq_ids == ["seq1", "seq2", "seq3"]
    consensus = cluster.consensus
    assert consensus is not None
    assert (
        str(consensus) == "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
    )


def test_clusters_class(sequences):
    cluster1 = Cluster("cluster0", sequences[:3])
    cluster2 = Cluster("cluster1", sequences[3:])
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


def test_cluster_mmseqs(fasta_file):
    cluster_info = cluster_mmseqs(
        fasta_file.name,
        threshold=0.99,
        temp_dir=tempfile.gettempdir(),
        as_dict=True,
    )
    assert len(cluster_info) == 1
    assert len(cluster_info["cluster0"]["seq_ids"]) == 3
    os.unlink(fasta_file.name)


def test_cluster_vsearch(fasta_file):
    cluster_info = cluster_vsearch(
        fasta_file.name,
        threshold=0.99,
        temp_dir=tempfile.gettempdir(),
        as_dict=True,
    )
    assert len(cluster_info) == 1
    assert len(cluster_info["cluster0"]["seq_ids"]) == 3
    os.unlink(fasta_file.name)


def test_cluster_with_list_of_sequences():
    seqs = ["ATCG", "ATCG", "CTAG", "CTAG"]
    clusters = cluster(seqs, threshold=0.9)
    assert len(clusters) == 2
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 2


def test_cluster_with_fasta_file():
    fasta_file = tempfile.NamedTemporaryFile(delete=False)
    fasta_file.write(b">seq1\nATCG\n>seq2\nATCG\n>seq3\nCTAG\n>seq4\nCTAG\n")
    fasta_file.close()
    clusters = cluster(fasta_file.name, threshold=0.9)
    assert len(clusters) == 2
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 2
    os.unlink(fasta_file.name)


def test_cluster_with_seqrecord_objects():
    seqs = [
        SeqRecord("ATCG"),
        SeqRecord("ATCG"),
        SeqRecord("CTAG"),
        SeqRecord("CTAG"),
    ]
    clusters = cluster(seqs, threshold=0.9)
    assert len(clusters) == 2
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 2


def test_cluster_with_sequence_objects():
    seqs = [Sequence("ATCG"), Sequence("ATCG"), Sequence("CTAG"), Sequence("CTAG")]
    clusters = cluster(seqs, threshold=0.9)
    assert len(clusters) == 2
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 2


def test_cluster_with_list_of_lists():
    seqs = [["seq1", "ATCG"], ["seq2", "ATCG"], ["seq3", "CTAG"], ["seq4", "CTAG"]]
    clusters = cluster(seqs, threshold=0.9)
    assert len(clusters) == 2
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 2
