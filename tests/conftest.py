"""Shared fixtures for the abutils test suite."""

import json
import os

import pytest

from abutils.core.lineage import Lineage
from abutils.core.pair import Pair
from abutils.core.sequence import Sequence

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


# ---------------------
#   Raw data loaders
# ---------------------


@pytest.fixture(scope="session")
def data_dir():
    """Path to the test data directory."""
    return DATA_DIR


@pytest.fixture(scope="session")
def pg9_json():
    """PG9 antibody annotation dict loaded from JSON."""
    with open(os.path.join(DATA_DIR, "PG9.json")) as f:
        return json.load(f)


@pytest.fixture(scope="session")
def pg9_fasta_path():
    """Path to the PG9 FASTA file."""
    return os.path.join(DATA_DIR, "PG9.fasta")


@pytest.fixture(scope="session")
def bcr_tsv_path():
    """Path to the BCR AIRR TSV file."""
    return os.path.join(DATA_DIR, "bcr_sequences.tsv")


@pytest.fixture(scope="session")
def tcr_tsv_path():
    """Path to the TCR AIRR TSV file."""
    return os.path.join(DATA_DIR, "tcr_sequences.tsv")


# ----------------------
#   Sequence fixtures
# ----------------------


@pytest.fixture(scope="session")
def pg9_sequence(pg9_json):
    """PG9 antibody as a Sequence object with full AIRR annotations."""
    return Sequence(pg9_json)


@pytest.fixture(scope="session")
def _bcr_sequence_dicts():
    """Load BCR sequence dicts from AIRR TSV."""
    import polars as pl

    df = pl.read_csv(os.path.join(DATA_DIR, "bcr_sequences.tsv"), separator="\t")
    return df.to_dicts()


@pytest.fixture(scope="session")
def _tcr_sequence_dicts():
    """Load TCR sequence dicts from AIRR TSV."""
    import polars as pl

    df = pl.read_csv(os.path.join(DATA_DIR, "tcr_sequences.tsv"), separator="\t")
    return df.to_dicts()


@pytest.fixture(scope="session")
def bcr_sequences(_bcr_sequence_dicts):
    """List of BCR Sequence objects from real antibody data."""
    return [Sequence(d) for d in _bcr_sequence_dicts]


@pytest.fixture(scope="session")
def tcr_sequences(_tcr_sequence_dicts):
    """List of TCR Sequence objects from real antibody data."""
    return [Sequence(d) for d in _tcr_sequence_dicts]


@pytest.fixture(scope="session")
def bcr_heavy_sequences(bcr_sequences):
    """BCR heavy chain sequences only."""
    return [s for s in bcr_sequences if s["locus"] == "IGH"]


@pytest.fixture(scope="session")
def bcr_light_sequences(bcr_sequences):
    """BCR light chain sequences only (kappa + lambda)."""
    return [s for s in bcr_sequences if s["locus"] in ("IGK", "IGL")]


@pytest.fixture(scope="session")
def tcr_alpha_sequences(tcr_sequences):
    """TCR alpha chain sequences only."""
    return [s for s in tcr_sequences if s["locus"] == "TRA"]


@pytest.fixture(scope="session")
def tcr_beta_sequences(tcr_sequences):
    """TCR beta chain sequences only."""
    return [s for s in tcr_sequences if s["locus"] == "TRB"]


# -------------------
#   Pair fixtures
# -------------------


@pytest.fixture(scope="session")
def bcr_pairs(bcr_heavy_sequences, bcr_light_sequences):
    """BCR Pair objects from real data, matched by index."""
    pairs = []
    for h, l in zip(bcr_heavy_sequences, bcr_light_sequences):
        pairs.append(Pair([h, l]))
    return pairs


@pytest.fixture(scope="session")
def tcr_pairs(tcr_alpha_sequences, tcr_beta_sequences):
    """TCR Pair objects from real data, matched by index."""
    pairs = []
    for a, b in zip(tcr_alpha_sequences, tcr_beta_sequences):
        pairs.append(Pair([a, b]))
    return pairs


# ----------------------
#   Lineage fixtures
# ----------------------


@pytest.fixture(scope="session")
def bcr_lineage_pairs(bcr_heavy_sequences, bcr_light_sequences):
    """BCR Pair objects annotated with lineage IDs for lineage grouping."""
    # Lineage 1: VRC01 (first 3 heavy + first 3 light)
    pairs = []
    for i in range(3):
        h = Sequence(
            {
                **bcr_heavy_sequences[i].annotations,
                "clonify": {"id": "lineage_001"},
            }
        )
        l = bcr_light_sequences[i]
        pairs.append(Pair([h, l]))
    # Lineage 2: PGT128 (last 2 heavy + last 2 light)
    for i in range(3, 5):
        h = Sequence(
            {
                **bcr_heavy_sequences[i].annotations,
                "clonify": {"id": "lineage_002"},
            }
        )
        l = bcr_light_sequences[i]
        pairs.append(Pair([h, l]))
    return pairs


@pytest.fixture(scope="session")
def bcr_lineage(bcr_lineage_pairs):
    """A Lineage object built from the first 3 BCR pairs (VRC01 lineage)."""
    return Lineage(bcr_lineage_pairs[:3])
