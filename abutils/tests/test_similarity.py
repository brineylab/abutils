import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from ..core.pair import Pair
from ..core.sequence import Sequence
from ..tools.similarity import (
    RepertoireSimilarity,
    RepertoireSimilarities,
    bray_curtis,
    cosine_similarity,
    jaccard_similarity,
    jensen_shannon,
    kullback_leibler,
    morisita_horn,
    renkonen,
    repertoire_similarity,
)


# ==================
#      FIXTURES
# ==================


@pytest.fixture
def sample_distributions():
    """Two sample frequency distributions for similarity testing."""
    sample1 = [10, 20, 30, 40, 0]
    sample2 = [15, 25, 25, 35, 5]
    return sample1, sample2


@pytest.fixture
def identical_distributions():
    """Identical distributions that should yield high similarity."""
    sample = [10, 20, 30, 40]
    return sample, sample.copy()


@pytest.fixture
def disjoint_distributions():
    """Disjoint distributions for testing edge cases."""
    sample1 = [10, 0, 20, 0]
    sample2 = [0, 15, 0, 25]
    return sample1, sample2


@pytest.fixture
def repertoire_sequences():
    """Two repertoires of Sequence objects for repertoire_similarity testing."""
    rep1 = [
        Sequence(
            {
                "sequence_id": f"rep1_seq{i}",
                "sequence": "ATCGATCG",
                "v_gene": f"IGHV{i % 3 + 1}-2",
                "j_gene": "IGHJ4",
                "cdr3_length": 12 + (i % 3),
                "locus": "IGH",
            }
        )
        for i in range(20)
    ]
    rep2 = [
        Sequence(
            {
                "sequence_id": f"rep2_seq{i}",
                "sequence": "ATCGATCG",
                "v_gene": f"IGHV{(i + 1) % 3 + 1}-2",
                "j_gene": "IGHJ4",
                "cdr3_length": 12 + ((i + 1) % 3),
                "locus": "IGH",
            }
        )
        for i in range(20)
    ]
    return rep1, rep2


@pytest.fixture
def similar_repertoire_sequences():
    """Two highly similar repertoires for testing."""
    rep1 = [
        Sequence(
            {
                "sequence_id": f"rep1_seq{i}",
                "sequence": "ATCGATCG",
                "v_gene": "IGHV1-2",
                "j_gene": "IGHJ4",
                "cdr3_length": 12,
                "locus": "IGH",
            }
        )
        for i in range(20)
    ]
    rep2 = [
        Sequence(
            {
                "sequence_id": f"rep2_seq{i}",
                "sequence": "ATCGATCG",
                "v_gene": "IGHV1-2",
                "j_gene": "IGHJ4",
                "cdr3_length": 12,
                "locus": "IGH",
            }
        )
        for i in range(20)
    ]
    return rep1, rep2


# ================================
#    INDIVIDUAL METRIC TESTS
# ================================


def test_morisita_horn_basic(sample_distributions):
    """Test morisita_horn returns a float between 0 and 1."""
    sample1, sample2 = sample_distributions
    result = morisita_horn(sample1, sample2)
    assert isinstance(result, float)
    assert 0 <= result <= 1


def test_morisita_horn_identical(identical_distributions):
    """Test morisita_horn returns ~1.0 for identical distributions."""
    sample1, sample2 = identical_distributions
    result = morisita_horn(sample1, sample2)
    assert result == pytest.approx(1.0, rel=1e-6)


def test_kullback_leibler_basic(sample_distributions):
    """Test kullback_leibler returns a float."""
    sample1, sample2 = sample_distributions
    # Filter out zeros for KL divergence
    s1 = [x for x in sample1 if x > 0]
    s2 = [x for x in sample2 if x > 0]
    # Need equal length arrays
    min_len = min(len(s1), len(s2))
    result = kullback_leibler(s1[:min_len], s2[:min_len])
    assert isinstance(result, float)


def test_jensen_shannon_basic(sample_distributions):
    """Test jensen_shannon returns a float between 0 and 1."""
    sample1, sample2 = sample_distributions
    result = jensen_shannon(sample1, sample2)
    assert isinstance(result, float)


def test_jaccard_similarity_basic(sample_distributions):
    """Test jaccard_similarity returns a float between 0 and 1."""
    sample1, sample2 = sample_distributions
    result = jaccard_similarity(sample1, sample2)
    assert isinstance(result, float)
    assert 0 <= result <= 1


def test_jaccard_similarity_identical(identical_distributions):
    """Test jaccard_similarity returns 1.0 for identical distributions."""
    sample1, sample2 = identical_distributions
    result = jaccard_similarity(sample1, sample2)
    assert result == pytest.approx(1.0, rel=1e-6)


def test_renkonen_basic(sample_distributions):
    """Test renkonen returns a float between 0 and 1."""
    sample1, sample2 = sample_distributions
    result = renkonen(sample1, sample2)
    assert isinstance(result, float)
    assert 0 <= result <= 1


def test_renkonen_identical(identical_distributions):
    """Test renkonen returns 1.0 for identical distributions."""
    sample1, sample2 = identical_distributions
    result = renkonen(sample1, sample2)
    assert result == pytest.approx(1.0, rel=1e-6)


def test_bray_curtis_basic(sample_distributions):
    """Test bray_curtis returns a float."""
    sample1, sample2 = sample_distributions
    result = bray_curtis(sample1, sample2)
    assert isinstance(result, float)


def test_bray_curtis_identical(identical_distributions):
    """Test bray_curtis returns 0.0 for identical distributions (distance metric)."""
    sample1, sample2 = identical_distributions
    result = bray_curtis(sample1, sample2)
    assert result == pytest.approx(0.0, abs=1e-6)


def test_cosine_similarity_basic(sample_distributions):
    """Test cosine_similarity returns a float between 0 and 1."""
    sample1, sample2 = sample_distributions
    result = cosine_similarity(sample1, sample2)
    assert isinstance(result, float)
    assert 0 <= result <= 1


def test_cosine_similarity_identical(identical_distributions):
    """Test cosine_similarity returns 1.0 for identical distributions."""
    sample1, sample2 = identical_distributions
    result = cosine_similarity(sample1, sample2)
    assert result == pytest.approx(1.0, rel=1e-6)


# ====================================
#    REPERTOIRE SIMILARITY CLASS
# ====================================


def test_repertoire_similarity_class_init():
    """Test RepertoireSimilarity initialization."""
    similarities = [0.8, 0.85, 0.82]
    rs = RepertoireSimilarity(
        similarities=similarities,
        method="morisita-horn",
        sample1_name="sample_a",
        sample2_name="sample_b",
    )
    assert rs.similarities == similarities
    assert rs.method == "morisita-horn"
    assert rs.sample1_name == "sample_a"
    assert rs.sample2_name == "sample_b"


def test_repertoire_similarity_default_names():
    """Test RepertoireSimilarity uses default names when not provided."""
    rs = RepertoireSimilarity(similarities=[0.8], method="jaccard")
    assert rs.sample1_name == "sample1"
    assert rs.sample2_name == "sample2"


def test_repertoire_similarity_len():
    """Test RepertoireSimilarity __len__ method."""
    similarities = [0.8, 0.85, 0.82, 0.79, 0.81]
    rs = RepertoireSimilarity(similarities=similarities, method="morisita-horn")
    assert len(rs) == 5


def test_repertoire_similarity_df_property():
    """Test RepertoireSimilarity df property returns a DataFrame."""
    similarities = [0.8, 0.85, 0.82]
    rs = RepertoireSimilarity(
        similarities=similarities,
        method="morisita-horn",
        sample1_name="sample_a",
        sample2_name="sample_b",
    )
    df = rs.df
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3
    assert "similarity" in df.columns
    assert "method" in df.columns
    assert "sample1" in df.columns
    assert "sample2" in df.columns
    assert "iteration" in df.columns


def test_repertoire_similarity_mean_median():
    """Test RepertoireSimilarity mean and median properties."""
    similarities = [0.8, 0.85, 0.82]
    rs = RepertoireSimilarity(similarities=similarities, method="morisita-horn")
    assert rs.mean == pytest.approx(np.mean(similarities), rel=1e-6)
    assert rs.median == pytest.approx(np.median(similarities), rel=1e-6)


def test_repertoire_similarity_std_sem():
    """Test RepertoireSimilarity std and sem properties."""
    similarities = [0.8, 0.85, 0.82, 0.79, 0.81]
    rs = RepertoireSimilarity(similarities=similarities, method="morisita-horn")
    assert rs.std == pytest.approx(np.std(similarities), rel=1e-6)
    expected_sem = np.std(similarities) / np.sqrt(len(similarities))
    assert rs.sem == pytest.approx(expected_sem, rel=1e-6)


def test_repertoire_similarity_min_max():
    """Test RepertoireSimilarity min and max properties."""
    similarities = [0.8, 0.85, 0.82]
    rs = RepertoireSimilarity(similarities=similarities, method="morisita-horn")
    assert rs.min == pytest.approx(0.8, rel=1e-6)
    assert rs.max == pytest.approx(0.85, rel=1e-6)


def test_repertoire_similarity_repr():
    """Test RepertoireSimilarity __repr__ method."""
    rs = RepertoireSimilarity(similarities=[0.8], method="morisita-horn")
    assert "morisita-horn" in repr(rs)


# =====================================
#    REPERTOIRE SIMILARITIES CLASS
# =====================================


def test_repertoire_similarities_init():
    """Test RepertoireSimilarities initialization."""
    rs1 = RepertoireSimilarity([0.8, 0.82], "morisita-horn", "a", "b")
    rs2 = RepertoireSimilarity([0.75, 0.78], "morisita-horn", "a", "c")
    rss = RepertoireSimilarities([rs1, rs2])
    assert len(rss) == 2


def test_repertoire_similarities_init_empty():
    """Test RepertoireSimilarities initialization with no arguments."""
    rss = RepertoireSimilarities()
    assert len(rss) == 0


def test_repertoire_similarities_add():
    """Test RepertoireSimilarities add method."""
    rss = RepertoireSimilarities()
    rs = RepertoireSimilarity([0.8], "morisita-horn", "a", "b")
    rss.add(rs)
    assert len(rss) == 1


def test_repertoire_similarities_add_operator():
    """Test RepertoireSimilarities __add__ method with RepertoireSimilarity."""
    rss = RepertoireSimilarities()
    rs = RepertoireSimilarity([0.8], "morisita-horn", "a", "b")
    rss = rss + rs
    assert len(rss) == 1


def test_repertoire_similarities_add_operator_with_similarities():
    """Test RepertoireSimilarities __add__ method with another RepertoireSimilarities."""
    rs1 = RepertoireSimilarity([0.8], "morisita-horn", "a", "b")
    rs2 = RepertoireSimilarity([0.75], "morisita-horn", "a", "c")
    rss1 = RepertoireSimilarities([rs1])
    rss2 = RepertoireSimilarities([rs2])
    rss_combined = rss1 + rss2
    assert len(rss_combined) == 2


def test_repertoire_similarities_iter():
    """Test RepertoireSimilarities __iter__ method."""
    rs1 = RepertoireSimilarity([0.8], "morisita-horn", "a", "b")
    rs2 = RepertoireSimilarity([0.75], "morisita-horn", "a", "c")
    rss = RepertoireSimilarities([rs1, rs2])
    items = list(rss)
    assert len(items) == 2
    assert rs1 in items
    assert rs2 in items


def test_repertoire_similarities_df_property():
    """Test RepertoireSimilarities df property returns a DataFrame."""
    rs1 = RepertoireSimilarity([0.8, 0.82], "morisita-horn", "a", "b")
    rs2 = RepertoireSimilarity([0.75, 0.78], "morisita-horn", "a", "c")
    rss = RepertoireSimilarities([rs1, rs2])
    df = rss.df
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4  # 2 similarities * 2 iterations each


def test_repertoire_similarities_means():
    """Test RepertoireSimilarities means property."""
    rs1 = RepertoireSimilarity([0.8, 0.82], "morisita-horn", "a", "b")
    rs2 = RepertoireSimilarity([0.75, 0.78], "morisita-horn", "a", "c")
    rss = RepertoireSimilarities([rs1, rs2])
    means = rss.means
    assert len(means) == 2
    assert means[0] == pytest.approx(np.mean([0.8, 0.82]), rel=1e-6)
    assert means[1] == pytest.approx(np.mean([0.75, 0.78]), rel=1e-6)


def test_repertoire_similarities_squareform():
    """Test RepertoireSimilarities squareform method."""
    rs1 = RepertoireSimilarity([0.8], "morisita-horn", "a", "b")
    rs2 = RepertoireSimilarity([0.75], "morisita-horn", "b", "a")
    rss = RepertoireSimilarities([rs1, rs2])
    sq = rss.squareform()
    assert isinstance(sq, pd.DataFrame)


def test_repertoire_similarities_repr():
    """Test RepertoireSimilarities __repr__ method."""
    rs1 = RepertoireSimilarity([0.8], "morisita-horn", "a", "b")
    rss = RepertoireSimilarities([rs1])
    assert "1 comparisons" in repr(rss)


# ========================================
#    REPERTOIRE_SIMILARITY FUNCTION
# ========================================


def test_repertoire_similarity_two_repertoires(repertoire_sequences):
    """Test repertoire_similarity with two repertoires returns a float."""
    rep1, rep2 = repertoire_sequences
    result = repertoire_similarity([rep1, rep2])
    assert isinstance(result, float)


def test_repertoire_similarity_multiple_iterations(repertoire_sequences):
    """Test repertoire_similarity with multiple iterations returns RepertoireSimilarity."""
    rep1, rep2 = repertoire_sequences
    result = repertoire_similarity([rep1, rep2], n_iters=5)
    assert isinstance(result, RepertoireSimilarity)
    assert len(result) == 5


def test_repertoire_similarity_with_names(repertoire_sequences):
    """Test repertoire_similarity with custom names."""
    rep1, rep2 = repertoire_sequences
    result = repertoire_similarity(
        [rep1, rep2], names=["repertoire_a", "repertoire_b"], n_iters=3
    )
    assert result.sample1_name == "repertoire_a"
    assert result.sample2_name == "repertoire_b"


def test_repertoire_similarity_different_methods(similar_repertoire_sequences):
    """Test repertoire_similarity with different similarity methods."""
    rep1, rep2 = similar_repertoire_sequences
    methods = [
        "morisita-horn",
        "jaccard",
        "bray-curtis",
        "renkonen",
        "cosine",
    ]
    for method in methods:
        result = repertoire_similarity([rep1, rep2], method=method)
        assert isinstance(result, float), f"Method {method} did not return a float"


def test_repertoire_similarity_three_repertoires(repertoire_sequences):
    """Test repertoire_similarity with three repertoires returns RepertoireSimilarities."""
    rep1, rep2 = repertoire_sequences
    # Create a third repertoire
    rep3 = [
        Sequence(
            {
                "sequence_id": f"rep3_seq{i}",
                "sequence": "ATCGATCG",
                "v_gene": f"IGHV{i % 2 + 1}-2",
                "j_gene": "IGHJ4",
                "cdr3_length": 13,
                "locus": "IGH",
            }
        )
        for i in range(20)
    ]
    result = repertoire_similarity([rep1, rep2, rep3])
    assert isinstance(result, RepertoireSimilarities)
    # With 3 repertoires, we get 9 pairwise comparisons (including self)
    assert len(result) == 9


def test_repertoire_similarity_custom_features(repertoire_sequences):
    """Test repertoire_similarity with custom features."""
    rep1, rep2 = repertoire_sequences
    result = repertoire_similarity([rep1, rep2], features=["v_gene"])
    assert isinstance(result, float)


def test_repertoire_similarity_subsample_size(repertoire_sequences):
    """Test repertoire_similarity with custom subsample size."""
    rep1, rep2 = repertoire_sequences
    result = repertoire_similarity([rep1, rep2], subsample_size=10)
    assert isinstance(result, float)


def test_repertoire_similarity_with_seed(repertoire_sequences):
    """Test repertoire_similarity produces reproducible results with seed."""
    rep1, rep2 = repertoire_sequences
    result1 = repertoire_similarity([rep1, rep2], n_iters=3, seed=42)
    result2 = repertoire_similarity([rep1, rep2], n_iters=3, seed=42)
    assert result1.similarities == result2.similarities
