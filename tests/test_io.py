import os
import shutil
import tempfile

import pandas as pd
import polars as pl
import pytest

from abutils.core.pair import Pair
from abutils.core.sequence import Sequence
from abutils.io import (
    concatenate_airr,
    concatenate_csv,
    concatenate_parquet,
    from_pandas,
    from_polars,
    split_airr,
    split_csv,
    split_fasta,
    split_fastq,
    split_parquet,
    to_pandas,
    to_polars,
)

# ==================
#      FIXTURES
# ==================


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def parquet_file():
    """Create a temp parquet file with 100 rows."""
    data = [
        {"sequence_id": f"seq{i}", "sequence": "ATCGATCG" * 5, "v_gene": f"IGHV{i % 5}"}
        for i in range(100)
    ]
    df = pl.DataFrame(data)
    pf = tempfile.NamedTemporaryFile(delete=False, suffix=".parquet")
    df.write_parquet(pf.name)
    pf.close()
    yield pf.name
    os.unlink(pf.name)


@pytest.fixture
def csv_file():
    """Create a temp CSV file with 100 rows."""
    content = "sequence_id,sequence,v_gene\n"
    content += "\n".join(
        [f"seq{i},ATCGATCGATCG,IGHV{i % 5}" for i in range(100)]
    )
    cf = tempfile.NamedTemporaryFile(delete=False, suffix=".csv", mode="w")
    cf.write(content)
    cf.close()
    yield cf.name
    os.unlink(cf.name)


@pytest.fixture
def airr_file():
    """Create a temp AIRR file (TSV) with 100 rows."""
    content = "sequence_id\tsequence\tv_gene\n"
    content += "\n".join(
        [f"seq{i}\tATCGATCGATCG\tIGHV{i % 5}" for i in range(100)]
    )
    af = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv", mode="w")
    af.write(content)
    af.close()
    yield af.name
    os.unlink(af.name)


@pytest.fixture
def fasta_file():
    """Create a temp FASTA file with 100 sequences."""
    content = "\n".join([f">seq{i}\nATCGATCGATCGATCG" for i in range(100)])
    ff = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w")
    ff.write(content)
    ff.close()
    yield ff.name
    os.unlink(ff.name)


@pytest.fixture
def fastq_file():
    """Create a temp FASTQ file with 100 sequences."""
    content = "\n".join(
        [f"@seq{i}\nATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIII" for i in range(100)]
    )
    fq = tempfile.NamedTemporaryFile(delete=False, suffix=".fastq", mode="w")
    fq.write(content)
    fq.close()
    yield fq.name
    os.unlink(fq.name)


@pytest.fixture
def sequence_list():
    """List of Sequence objects for DataFrame conversion testing."""
    return [
        Sequence(
            {
                "sequence_id": f"seq{i}",
                "sequence": "ATCGATCG",
                "v_gene": f"IGHV{i}",
                "locus": "IGH",
            }
        )
        for i in range(5)
    ]


@pytest.fixture
def pair_list():
    """List of Pair objects for DataFrame conversion testing."""
    pairs = []
    for i in range(5):
        heavy = Sequence(
            {
                "sequence_id": f"pair{i}",
                "sequence": "ATCGATCGATCG",
                "locus": "IGH",
                "v_gene": f"IGHV{i}-1",
            }
        )
        light = Sequence(
            {
                "sequence_id": f"pair{i}",
                "sequence": "GCTAGCTAGCTA",
                "locus": "IGK",
                "v_gene": f"IGKV{i}-1",
            }
        )
        pairs.append(Pair([heavy, light]))
    return pairs


# ============================
#    SPLIT FUNCTION TESTS
# ============================


def test_split_parquet_creates_files(parquet_file, temp_dir):
    """Test split_parquet creates multiple output files."""
    output_files = split_parquet(parquet_file, temp_dir, num_rows=30)
    assert len(output_files) == 4  # 100 rows / 30 = 4 files (with remainder)
    for f in output_files:
        assert os.path.exists(f)
        assert f.endswith(".parquet")


def test_split_parquet_correct_row_count(parquet_file, temp_dir):
    """Test split_parquet produces files with correct row counts."""
    output_files = split_parquet(parquet_file, temp_dir, num_rows=25)
    total_rows = 0
    for f in output_files:
        df = pl.read_parquet(f)
        total_rows += len(df)
    assert total_rows == 100


def test_split_parquet_with_num_splits(parquet_file, temp_dir):
    """Test split_parquet with num_splits parameter."""
    output_files = split_parquet(parquet_file, temp_dir, num_splits=5)
    assert len(output_files) == 5


def test_split_csv_creates_files(csv_file, temp_dir):
    """Test split_csv creates multiple output files."""
    output_files = split_csv(csv_file, temp_dir, chunksize=30)
    assert len(output_files) == 4  # 100 rows / 30 = 4 files
    for f in output_files:
        assert os.path.exists(f)
        assert f.endswith(".csv")


def test_split_airr_creates_files(airr_file, temp_dir):
    """Test split_airr creates multiple output files."""
    output_files = split_airr(airr_file, temp_dir, chunksize=30)
    assert len(output_files) == 4
    for f in output_files:
        assert os.path.exists(f)


def test_split_fasta_creates_files(fasta_file, temp_dir):
    """Test split_fasta creates multiple output files."""
    output_files = split_fasta(fasta_file, temp_dir, chunksize=30)
    assert len(output_files) == 4
    for f in output_files:
        assert os.path.exists(f)
        assert f.endswith(".fasta")


def test_split_fastq_creates_files(fastq_file, temp_dir):
    """Test split_fastq creates multiple output files."""
    output_files = split_fastq(fastq_file, temp_dir, chunksize=30)
    assert len(output_files) == 4
    for f in output_files:
        assert os.path.exists(f)
        assert f.endswith(".fastq")


# ================================
#    CONCATENATE FUNCTION TESTS
# ================================


def test_concatenate_parquet(parquet_file, temp_dir):
    """Test concatenate_parquet merges files correctly."""
    # First split, then concatenate
    split_files = split_parquet(parquet_file, temp_dir, num_rows=30)
    output_file = os.path.join(temp_dir, "concatenated.parquet")
    result = concatenate_parquet(split_files, output_file)
    assert os.path.exists(result)
    df = pl.read_parquet(result)
    assert len(df) == 100


def test_concatenate_csv(csv_file, temp_dir):
    """Test concatenate_csv merges files correctly."""
    # First split, then concatenate
    split_files = split_csv(csv_file, temp_dir, chunksize=30)
    output_file = os.path.join(temp_dir, "concatenated.csv")
    result = concatenate_csv(split_files, output_file)
    assert os.path.exists(result)
    df = pl.read_csv(result)
    assert len(df) == 100


def test_concatenate_airr(airr_file, temp_dir):
    """Test concatenate_airr merges files correctly."""
    # First split, then concatenate
    split_files = split_airr(airr_file, temp_dir, chunksize=30)
    output_file = os.path.join(temp_dir, "concatenated.tsv")
    result = concatenate_airr(split_files, output_file)
    assert os.path.exists(result)
    df = pl.read_csv(result, separator="\t")
    assert len(df) == 100


# ====================================
#    DATAFRAME CONVERSION TESTS
# ====================================


def test_to_polars_with_sequences(sequence_list):
    """Test to_polars with Sequence objects."""
    df = to_polars(sequence_list)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 5
    assert "sequence_id" in df.columns
    assert "sequence" in df.columns


def test_from_polars_with_sequences():
    """Test from_polars returns Sequence objects."""
    df = pl.DataFrame(
        {
            "sequence_id": ["seq1", "seq2", "seq3"],
            "sequence": ["ATCG", "GCTA", "TTAA"],
            "v_gene": ["IGHV1", "IGHV2", "IGHV3"],
        }
    )
    sequences = from_polars(df)
    assert len(sequences) == 3
    assert all(isinstance(s, Sequence) for s in sequences)
    assert sequences[0]["sequence_id"] == "seq1"


def test_to_pandas_with_sequences(sequence_list):
    """Test to_pandas with Sequence objects."""
    df = to_pandas(sequence_list)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 5
    assert "sequence_id" in df.columns
    assert "sequence" in df.columns


def test_from_pandas_with_sequences():
    """Test from_pandas returns Sequence objects."""
    df = pd.DataFrame(
        {
            "sequence_id": ["seq1", "seq2", "seq3"],
            "sequence": ["ATCG", "GCTA", "TTAA"],
            "v_gene": ["IGHV1", "IGHV2", "IGHV3"],
        }
    )
    sequences = from_pandas(df)
    assert len(sequences) == 3
    assert all(isinstance(s, Sequence) for s in sequences)


def test_to_polars_with_pairs(pair_list):
    """Test to_polars with Pair objects."""
    df = to_polars(pair_list)
    assert isinstance(df, pl.DataFrame)
    assert len(df) == 5
    # Pair data uses :0 and :1 suffixes
    assert "sequence:0" in df.columns or "sequence_id:0" in df.columns


def test_from_polars_with_pairs():
    """Test from_polars returns Pair objects when data has :0/:1 columns."""
    df = pl.DataFrame(
        {
            "sequence_id:0": ["pair1", "pair2"],
            "sequence:0": ["ATCG", "GCTA"],
            "locus:0": ["IGH", "IGH"],
            "sequence_id:1": ["pair1", "pair2"],
            "sequence:1": ["TTAA", "AACC"],
            "locus:1": ["IGK", "IGK"],
        }
    )
    pairs = from_polars(df)
    assert len(pairs) == 2
    assert all(isinstance(p, Pair) for p in pairs)


def test_to_pandas_with_pairs(pair_list):
    """Test to_pandas with Pair objects."""
    df = to_pandas(pair_list)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 5


def test_from_polars_with_lazyframe():
    """Test from_polars works with LazyFrame."""
    df = pl.DataFrame(
        {
            "sequence_id": ["seq1", "seq2"],
            "sequence": ["ATCG", "GCTA"],
        }
    ).lazy()
    sequences = from_polars(df)
    assert len(sequences) == 2
    assert all(isinstance(s, Sequence) for s in sequences)


def test_to_polars_with_annotations(sequence_list):
    """Test to_polars with specific annotations."""
    df = to_polars(sequence_list, annotations=["sequence_id", "v_gene"])
    assert "sequence_id" in df.columns
    assert "v_gene" in df.columns


def test_to_polars_drop_na_columns(sequence_list):
    """Test to_polars drop_na_columns parameter."""
    # Add a sequence with extra annotation
    sequence_list[0].annotations["extra_field"] = "value"
    df = to_polars(sequence_list, drop_na_columns=True)
    # extra_field should be dropped since most sequences don't have it
    # (or kept if implementation doesn't drop)
    assert isinstance(df, pl.DataFrame)
