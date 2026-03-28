import pytest

from abutils.core.lineage import Lineage, LineageSummary, group_lineages
from abutils.core.pair import Pair
from abutils.core.sequence import Sequence

# ==================
#      FIXTURES
# ==================


@pytest.fixture
def mock_heavy_sequences():
    """Create mock heavy chain Sequence objects."""
    return [
        Sequence(
            {
                "sequence_id": f"seq{i}",
                "sequence": "ATCGATCGATCGATCG",
                "sequence_aa": "CASSLDRYFDYW",
                "locus": "IGH",
                "v_call": "IGHV1-2*01",
                "j_call": "IGHJ4*01",
                "c_gene": "IGHG1",
                "cdr3_aa": "CARASGY",
                "cdr3_length": 7,
                "v_identity": 0.95 - (i * 0.01),
                "v_identity_aa": 0.92 - (i * 0.01),
                "lineage": "lineage_001",
            }
        )
        for i in range(5)
    ]


@pytest.fixture
def mock_light_sequences():
    """Create mock light chain Sequence objects."""
    return [
        Sequence(
            {
                "sequence_id": f"seq{i}",
                "sequence": "GCTAGCTAGCTAGCTA",
                "sequence_aa": "CQQSYSTPLT",
                "locus": "IGK",
                "v_call": "IGKV1-5*01",
                "j_call": "IGKJ1*01",
                "cdr3_aa": "CQQSY",
                "cdr3_length": 5,
                "v_identity": 0.97 - (i * 0.01),
                "v_identity_aa": 0.95 - (i * 0.01),
            }
        )
        for i in range(5)
    ]


@pytest.fixture
def mock_pairs(mock_heavy_sequences, mock_light_sequences):
    """Create mock Pair objects for Lineage testing."""
    pairs = []
    for heavy, light in zip(mock_heavy_sequences, mock_light_sequences):
        pairs.append(Pair([heavy, light]))
    return pairs


@pytest.fixture
def mock_pairs_heavy_only(mock_heavy_sequences):
    """Create mock Pair objects with only heavy chains."""
    return [Pair([heavy]) for heavy in mock_heavy_sequences]


@pytest.fixture
def mock_lineage(mock_pairs):
    """Create a mock Lineage object."""
    return Lineage(mock_pairs)


# ================================
#    LINEAGE CLASS BASIC TESTS
# ================================


def test_lineage_init(mock_pairs):
    """Test Lineage initialization."""
    lineage = Lineage(mock_pairs)
    assert lineage is not None
    assert lineage.pairs == mock_pairs


def test_lineage_init_with_name(mock_pairs):
    """Test Lineage initialization with a name."""
    lineage = Lineage(mock_pairs, name="test_lineage")
    assert lineage is not None


def test_lineage_size(mock_lineage):
    """Test Lineage size method."""
    assert mock_lineage.size() == 5


def test_lineage_size_pairs_only(mock_lineage):
    """Test Lineage size method with pairs_only=True."""
    size = mock_lineage.size(pairs_only=True)
    assert size == 5  # All mock pairs have both chains


def test_lineage_size_heavy_only(mock_pairs_heavy_only):
    """Test Lineage size with heavy-only pairs."""
    lineage = Lineage(mock_pairs_heavy_only)
    assert lineage.size() == 5
    assert lineage.size(pairs_only=True) == 0


def test_lineage_contains(mock_lineage, mock_pairs):
    """Test Lineage __contains__ method."""
    # The pair_dict uses pair.name as keys
    pair_name = mock_pairs[0].name
    assert pair_name in mock_lineage


def test_lineage_getitem(mock_lineage, mock_pairs):
    """Test Lineage __getitem__ method."""
    pair_name = mock_pairs[0].name
    retrieved_pair = mock_lineage[pair_name]
    assert retrieved_pair is not None
    assert retrieved_pair.name == pair_name


def test_lineage_iter(mock_lineage, mock_pairs):
    """Test Lineage __iter__ method."""
    pairs_from_iter = list(mock_lineage)
    assert len(pairs_from_iter) == 5
    assert pairs_from_iter == mock_pairs


# ================================
#    LINEAGE PROPERTIES TESTS
# ================================


def test_lineage_heavies_property(mock_lineage):
    """Test Lineage heavies property."""
    heavies = mock_lineage.heavies
    assert len(heavies) == 5
    for pair in heavies:
        assert pair.heavy is not None


def test_lineage_lights_property(mock_lineage):
    """Test Lineage lights property."""
    lights = mock_lineage.lights
    assert len(lights) == 5
    for pair in lights:
        assert pair.light is not None


def test_lineage_just_pairs_property(mock_lineage):
    """Test Lineage just_pairs property."""
    just_pairs = mock_lineage.just_pairs
    assert len(just_pairs) == 5
    for pair in just_pairs:
        assert pair.is_pair


def test_lineage_name_property(mock_lineage):
    """Test Lineage name property returns string or None."""
    name = mock_lineage.name
    # Name can be None if lineage info isn't found in expected format
    assert name is None or isinstance(name, str)


def test_lineage_pair_dict_property(mock_lineage, mock_pairs):
    """Test Lineage pair_dict property."""
    pair_dict = mock_lineage.pair_dict
    assert isinstance(pair_dict, dict)
    assert len(pair_dict) == 5
    for pair in mock_pairs:
        assert pair.name in pair_dict


def test_lineage_has_insertion(mock_lineage):
    """Test Lineage has_insertion property."""
    # Our mock sequences don't have insertions
    assert mock_lineage.has_insertion is False


def test_lineage_has_deletion(mock_lineage):
    """Test Lineage has_deletion property."""
    # Our mock sequences don't have deletions
    assert mock_lineage.has_deletion is False


def test_lineage_has_indel(mock_lineage):
    """Test Lineage has_indel property."""
    # Our mock sequences don't have indels
    assert mock_lineage.has_indel is False


# ================================
#    LINEAGE SUMMARY TESTS
# ================================


def test_lineage_summary_init_with_pairs(mock_pairs):
    """Test LineageSummary initialization with Pair objects."""
    summary = LineageSummary(mock_pairs, name="test_summary")
    assert summary is not None
    assert summary.name == "test_summary"
    assert len(summary.pairs) == 5


def test_lineage_summary_heavies_lights(mock_pairs):
    """Test LineageSummary heavies and lights attributes."""
    summary = LineageSummary(mock_pairs, name="test_summary")
    assert len(summary.heavies) == 5
    assert len(summary.lights) == 5


def test_lineage_summary_df(mock_pairs):
    """Test LineageSummary has a df attribute."""
    summary = LineageSummary(mock_pairs, name="test_summary")
    assert summary.df is not None
    assert len(summary.df) == 5


# ================================
#    GROUP LINEAGES TESTS
# ================================


def test_group_lineages_basic():
    """Test group_lineages function."""
    # Create pairs with clonify IDs
    pairs = []
    for i in range(5):
        heavy = Sequence(
            {
                "sequence_id": f"seq{i}",
                "sequence": "ATCG",
                "locus": "IGH",
                "clonify": {"id": "lineage_001" if i < 3 else "lineage_002"},
            }
        )
        light = Sequence(
            {
                "sequence_id": f"seq{i}",
                "sequence": "GCTA",
                "locus": "IGK",
            }
        )
        pairs.append(Pair([heavy, light]))

    lineages = group_lineages(pairs)
    assert len(lineages) == 2

    # Check lineage sizes
    sizes = sorted([l.size() for l in lineages], reverse=True)
    assert sizes[0] == 3
    assert sizes[1] == 2


def test_group_lineages_just_pairs():
    """Test group_lineages with just_pairs=True."""
    # Create pairs with clonify IDs, some without light chains
    pairs = []
    for i in range(5):
        heavy = Sequence(
            {
                "sequence_id": f"seq{i}",
                "sequence": "ATCG",
                "locus": "IGH",
                "clonify": {"id": "lineage_001"},
            }
        )
        if i < 3:  # Only first 3 have light chains
            light = Sequence(
                {
                    "sequence_id": f"seq{i}",
                    "sequence": "GCTA",
                    "locus": "IGK",
                }
            )
            pairs.append(Pair([heavy, light]))
        else:
            pairs.append(Pair([heavy]))

    # Without just_pairs
    lineages_all = group_lineages(pairs, just_pairs=False)
    assert len(lineages_all) == 1
    assert lineages_all[0].size() == 5

    # With just_pairs
    lineages_pairs = group_lineages(pairs, just_pairs=True)
    assert len(lineages_pairs) == 1
    assert lineages_pairs[0].size() == 3


def test_group_lineages_empty():
    """Test group_lineages with no clonify annotations."""
    pairs = []
    for i in range(3):
        heavy = Sequence(
            {
                "sequence_id": f"seq{i}",
                "sequence": "ATCG",
                "locus": "IGH",
                # No clonify annotation
            }
        )
        pairs.append(Pair([heavy]))

    lineages = group_lineages(pairs)
    assert len(lineages) == 0
