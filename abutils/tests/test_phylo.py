import os
import platform
import pytest
import tempfile

import numpy as np

from matplotlib.path import Path

from ..core.sequence import Sequence
from ..tools.phylo import fasttree, phylogeny, Phylogeny, align_marker


@pytest.fixture
def sequences():
    seq1 = Sequence(
        {
            "sequence": "ATCGATCGATCGATCGGGCGATCGATCG",
            "seq": "ATCGATCGATCGATCGGGCGATCGATCG",
            "sequence_id": "seq1",
            "name": "seq1",
        }
    )
    seq2 = Sequence(
        {
            "sequence": "ATCGATCGATCGATCGGGGGATCGATCG",
            "seq": "ATCGATCGATCGATCGGGGGATCGATCG",
            "sequence_id": "seq2",
            "name": "seq2",
        }
    )
    seq3 = Sequence(
        {
            "sequence": "ATCGCATGATCGTTCGGGCGATCGATCG",
            "seq": "ATCGCATGATCGTTCGGGCGATCGATCG",
            "sequence_id": "seq3",
            "name": "seq3",
        }
    )
    return [seq1, seq2, seq3]


@pytest.fixture
def fasta_string(sequences):
    return "\n".join([s.fasta for s in sequences])


@pytest.fixture
def fasta_file(fasta_string):
    fasta_file = tempfile.NamedTemporaryFile(delete=False)
    fasta_file.write(fasta_string.encode("utf-8"))
    fasta_file.close()
    return fasta_file


@pytest.fixture
def sequences_aa():
    seq1 = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVISPLPSQAMDDLMLSPDDIEQ", id="seq1")
    seq2 = Sequence("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQ", id="seq2")
    seq3 = Sequence("MEEPQSDPSVEPPLSQETFSDIWKLLPENNVISPLPSQAMDDLMLSPDDIEQ", id="seq3")
    return [seq1, seq2, seq3]


@pytest.fixture
def fasta_string_aa(sequences_aa):
    return "\n".join([s.fasta for s in sequences_aa])


@pytest.fixture
def fasta_file_aa(fasta_string_aa):
    fasta_file = tempfile.NamedTemporaryFile(delete=False)
    fasta_file.write(fasta_string_aa.encode("utf-8"))
    fasta_file.close()
    return fasta_file


@pytest.fixture
def root_sequence():
    return Sequence(seq="ATCGATCGATCGATCGGGCGATCGATCG", id="root")


# ---------------------------
#       fasttree
# ---------------------------``


def test_fasttree_with_file_input(fasta_file):
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    tree_path = fasttree(
        fasta_file.name, tree_file=tree_file.name, debug=True, quiet=False
    )
    assert os.path.isfile(tree_path)
    with open(tree_path, "r") as f:
        tree_string = f.read()
    assert tree_string.startswith("(")
    os.remove(tree_path)


def test_fasttree_with_string_input(fasta_string):
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    tree_path = fasttree(
        fasta_string, tree_file=tree_file.name, debug=True, quiet=False
    )
    assert os.path.isfile(tree_path)
    with open(tree_path, "r") as f:
        tree_string = f.read()
    assert tree_string.startswith("(")
    os.remove(tree_path)


def test_fasttree_with_aa_file_input(fasta_file_aa):
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    tree_path = fasttree(
        fasta_file_aa.name,
        tree_file=tree_file.name,
        is_aa=True,
        debug=True,
        quiet=False,
    )
    assert os.path.isfile(tree_path)
    with open(tree_path, "r") as f:
        tree_string = f.read()
    assert tree_string.startswith("(")
    os.remove(tree_path)


def test_fasttree_with_aa_string_input(fasta_string_aa):
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    tree_path = fasttree(
        fasta_string_aa, tree_file=tree_file.name, is_aa=True, debug=True, quiet=False
    )
    assert os.path.isfile(tree_path)
    with open(tree_path, "r") as f:
        tree_string = f.read()
    assert tree_string.startswith("(")
    os.remove(tree_path)


def test_fasttree_with_custom_binary(fasta_file):
    test_bin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fasttree_bin = os.path.join(
        test_bin_dir, f"tests/bin/fasttree_{platform.system().lower()}"
    )
    tree_file = tempfile.NamedTemporaryFile(delete=False)
    tree_path = fasttree(
        fasta_file.name,
        tree_file=tree_file.name,
        fasttree_bin=fasttree_bin,
        debug=True,
        quiet=False,
    )
    assert os.path.isfile(tree_path)
    with open(tree_path, "r") as f:
        tree_string = f.read()
    assert tree_string.startswith("(")
    os.remove(tree_path)


# ---------------------------
#        phylogeny
# ---------------------------


def test_phylogeny_with_sequence_input(sequences):
    phylo = phylogeny(sequences)
    assert isinstance(phylo, Phylogeny)
    assert len(phylo.sequences) == 3


def test_phylogeny_with_file_input(fasta_file):
    phylo = phylogeny(fasta_file.name)
    assert isinstance(phylo, Phylogeny)
    assert len(phylo.sequences) == 3


def test_phylogeny_with_list_input(sequences):
    phylo = phylogeny(sequences)
    assert isinstance(phylo, Phylogeny)
    assert len(phylo.sequences) == 3


def test_phylogeny_with_name(fasta_file):
    name = "test_phylogeny_with_name"
    phylo = phylogeny(fasta_file.name, name=name)
    assert isinstance(phylo, Phylogeny)
    assert phylo.name == name


def test_phylogeny_with_root_sequence(fasta_file, root_sequence):
    phylo = phylogeny(fasta_file.name, root=root_sequence)
    assert isinstance(phylo, Phylogeny)
    assert phylo.root.id == root_sequence.id


def test_phylogeny_with_root_sequence_id(fasta_file):
    root_seq_id = "seq1"
    phylo = phylogeny(fasta_file.name, root=root_seq_id)
    assert isinstance(phylo, Phylogeny)
    assert phylo.root.id == root_seq_id


def test_phylogeny_with_cluster(fasta_file):
    phylo = phylogeny(fasta_file.name, cluster=True, clustering_algo="cdhit")
    assert isinstance(phylo, Phylogeny)
    assert len(phylo.sequences) == 3
    assert len(phylo.clusters) == 3


def test_phylogeny_with_clustering_threshold(fasta_file):
    clustering_threshold = 0.95
    phylo = phylogeny(
        fasta_file.name,
        cluster=True,
        clustering_threshold=clustering_threshold,
        clustering_algo="cdhit",
    )
    assert isinstance(phylo, Phylogeny)
    assert len(phylo.sequences) == 3
    assert len(phylo.clusters) == 2


def test_phylogeny_with_rename_dict(fasta_file):
    rename_dict = {"seq1": "new_seq1", "seq2": "new_seq2"}
    phylo = phylogeny(fasta_file.name, rename=rename_dict)
    assert isinstance(phylo, Phylogeny)
    phylo_ids = [s.id for s in phylo.sequences]
    assert "new_seq1" in phylo_ids
    assert "new_seq2" in phylo_ids
    assert "seq1" not in phylo_ids
    assert "seq2" not in phylo_ids


def test_phylogeny_with_rename_callable(fasta_file):
    def rename_func(old_name):
        if old_name == "seq1":
            return "new_seq1"
        else:
            return None

    phylo = phylogeny(fasta_file.name, rename=rename_func)
    assert isinstance(phylo, Phylogeny)
    phylo_ids = [s.id for s in phylo.sequences]
    assert "new_seq1" in phylo_ids
    assert "seq1" not in phylo_ids


def test_phylogeny_with_id_key(sequences):
    phylo = phylogeny(sequences, id_key="name")
    assert isinstance(phylo, Phylogeny)
    assert len(phylo.sequences) == 3


def test_phylogeny_with_sequence_key(sequences):
    phylo = phylogeny(sequences, sequence_key="seq")
    assert isinstance(phylo, Phylogeny)
    assert len(phylo.sequences) == 3


# ---------------------------
#      align_marker
# ---------------------------


def test_align_marker_with_center_alignment():
    marker = "o"
    halign = "center"
    valign = "middle"
    marker_path = align_marker(marker, halign=halign, valign=valign)
    assert isinstance(marker_path, Path)
    assert np.allclose(marker_path.vertices[:, 0].max(), 0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 0].min(), -0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].max(), 0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].min(), -0.5, atol=0.05)
    assert np.allclose(marker_path.vertices.mean(axis=0), [0, 0], atol=0.05)


def test_align_marker_with_left_alignment():
    marker = "o"
    halign = "left"
    valign = "middle"
    marker_path = align_marker(marker, halign=halign, valign=valign)
    assert isinstance(marker_path, Path)
    assert np.allclose(marker_path.vertices[:, 0].max(), 1, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 0].min(), 0, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].max(), 0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].min(), -0.5, atol=0.05)
    assert np.allclose(marker_path.vertices.mean(axis=0), [0.5, 0], atol=0.05)


def test_align_marker_with_right_alignment():
    marker = "o"
    halign = "right"
    valign = "middle"
    marker_path = align_marker(marker, halign=halign, valign=valign)
    assert isinstance(marker_path, Path)
    assert np.allclose(marker_path.vertices[:, 0].max(), 0, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 0].min(), -1, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].max(), 0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].min(), -0.5, atol=0.05)
    assert np.allclose(marker_path.vertices.mean(axis=0), [-0.5, 0], atol=0.05)


def test_align_marker_with_top_alignment():
    marker = "o"
    halign = "center"
    valign = "top"
    marker_path = align_marker(marker, halign=halign, valign=valign)
    assert isinstance(marker_path, Path)
    assert np.allclose(marker_path.vertices[:, 0].max(), 0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 0].min(), -0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].max(), 0, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].min(), -1, atol=0.05)
    assert np.allclose(marker_path.vertices.mean(axis=0), [0, -0.5], atol=0.05)


def test_align_marker_with_bottom_alignment():
    marker = "o"
    halign = "center"
    valign = "bottom"
    marker_path = align_marker(marker, halign=halign, valign=valign)
    assert isinstance(marker_path, Path)
    assert np.allclose(marker_path.vertices[:, 0].max(), 0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 0].min(), -0.5, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].max(), 1, atol=0.05)
    assert np.allclose(marker_path.vertices[:, 1].min(), 0, atol=0.05)
    assert np.allclose(marker_path.vertices.mean(axis=0), [0, 0.5], atol=0.05)
