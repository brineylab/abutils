#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import tempfile
from pathlib import Path
from unittest import mock

import polars as pl
import pytest

from abutils import Sequence
from abutils.tools import preprocessing


# Fixtures
@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    dir_path = tempfile.mkdtemp(prefix="tempdir")
    yield dir_path
    shutil.rmtree(dir_path)


@pytest.fixture
def fastq_file(temp_dir):
    """Create a simple FASTQ file for testing."""
    fastq_path = os.path.join(temp_dir, "test_file.fastq")
    with open(fastq_path, "w") as f:
        f.write("test content")
    return fastq_path


@pytest.fixture
def illumina_fastq(temp_dir):
    """Create a FASTQ file with Illumina naming convention."""
    fastq_path = os.path.join(temp_dir, "Sample1_S1_L001_R1_001.fastq")
    with open(fastq_path, "w") as f:
        f.write("@header1\nACGT\n+\nIIII\n")
    return fastq_path


@pytest.fixture
def element_fastq(temp_dir):
    """Create a FASTQ file with Element naming convention."""
    fastq_path = os.path.join(temp_dir, "Sample1_R1.fastq")
    with open(fastq_path, "w") as f:
        f.write("@header1\nACGT\n+\nIIII\n")
    return fastq_path


@pytest.fixture
def illumina_pair(temp_dir):
    """Create a pair of Illumina FASTQ files."""
    r1_path = os.path.join(temp_dir, "Sample1_S1_L001_R1_001.fastq")
    r2_path = os.path.join(temp_dir, "Sample1_S1_L001_R2_001.fastq")

    with open(r1_path, "w") as f:
        f.write("@header1\nACGT\n+\nIIII\n")

    with open(r2_path, "w") as f:
        f.write("@header1\nTGCA\n+\nIIII\n")

    return r1_path, r2_path


@pytest.fixture
def merged_dir(temp_dir):
    """Create an output directory for merged files."""
    merged_dir = os.path.join(temp_dir, "merged")
    os.makedirs(merged_dir, exist_ok=True)
    return merged_dir


@pytest.fixture
def log_dir(temp_dir):
    """Create a directory for log files."""
    log_dir = os.path.join(temp_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    return log_dir


@pytest.fixture
def multi_sample_fastqs(temp_dir):
    """Create multiple sample FASTQ files."""
    r1_path = os.path.join(temp_dir, "Sample1_S1_L001_R1_001.fastq")
    r2_path = os.path.join(temp_dir, "Sample1_S1_L001_R2_001.fastq")
    r1_path2 = os.path.join(temp_dir, "Sample2_S2_L001_R1_001.fastq")
    r2_path2 = os.path.join(temp_dir, "Sample2_S2_L001_R2_001.fastq")

    for path in [r1_path, r2_path, r1_path2, r2_path2]:
        with open(path, "w") as f:
            f.write("@header1\nACGT\n+\nIIII\n")

    return {
        "dirname": temp_dir,
        "r1_path": r1_path,
        "r2_path": r2_path,
        "r1_path2": r1_path2,
        "r2_path2": r2_path2,
    }


@pytest.fixture
def airr_files(temp_dir):
    """Create temporary AIRR files for testing."""
    sample1_path = os.path.join(temp_dir, "sample1.tsv")
    sample2_path = os.path.join(temp_dir, "sample2.tsv")

    # Create a simple AIRR dataframe
    data1 = {
        "sequence_id": ["seq1", "seq2", "seq3"],
        "sequence": ["ACGT", "ACGT", "TGCA"],
        "umi": ["UMI1", "UMI2", "UMI3"],
        "v_gene": ["IGHV1-1*01", "IGHV1-2*01", "IGHV2-1*01"],
        "j_gene": ["IGHJ1*01", "IGHJ2*01", "IGHJ1*01"],
        "locus": ["IGH", "IGH", "IGH"],
    }

    data2 = {
        "sequence_id": ["seq4", "seq5", "seq6"],
        "sequence": ["ACGT", "GTCA", "TGCA"],
        "umi": ["UMI4", "UMI5", "UMI6"],
        "v_gene": ["IGHV1-1*01", "IGHV1-2*01", "IGHV2-1*01"],
        "j_gene": ["IGHJ1*01", "IGHJ2*01", "IGHJ1*01"],
        "locus": ["IGH", "IGH", "IGH"],
    }

    df1 = pl.DataFrame(data1)
    df2 = pl.DataFrame(data2)

    df1.write_csv(sample1_path, separator="\t")
    df2.write_csv(sample2_path, separator="\t")

    # Create output directory
    output_dir = os.path.join(temp_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    return {
        "dirname": temp_dir,
        "sample1_path": sample1_path,
        "sample2_path": sample2_path,
        "output_dir": output_dir,
    }


# Tests for FASTQFile
def test_fastq_file_init(fastq_file):
    """Test initialization of FASTQFile."""
    fastq_obj = preprocessing.FASTQFile(fastq_file)
    assert fastq_obj.file == fastq_file
    assert fastq_obj.path == os.path.abspath(fastq_file)
    assert fastq_obj.basename == "test_file.fastq"
    assert fastq_obj.dir == os.path.dirname(os.path.abspath(fastq_file))
    assert fastq_obj.filename == "test_file"


# Tests for IlluminaFile
def test_illumina_file_init(illumina_fastq):
    """Test initialization of IlluminaFile."""
    illumina_file = preprocessing.IlluminaFile(illumina_fastq)
    assert illumina_file.schema == "illumina"
    assert illumina_file.filename == "Sample1_S1_L001_R1_001"


def test_illumina_file_properties(illumina_fastq):
    """Test properties of IlluminaFile."""
    illumina_file = preprocessing.IlluminaFile(illumina_fastq)
    assert illumina_file.name == "Sample1"
    assert illumina_file.sample == "S1"
    assert illumina_file.lane == "L001"
    assert illumina_file.read == "R1"
    assert illumina_file.number == "001"


# Tests for ElementFile
def test_element_file_init(element_fastq):
    """Test initialization of ElementFile."""
    element_file = preprocessing.ElementFile(element_fastq)
    assert element_file.schema == "element"
    assert element_file.filename == "Sample1_R1"


def test_element_file_properties(element_fastq):
    """Test properties of ElementFile."""
    element_file = preprocessing.ElementFile(element_fastq)
    assert element_file.name == "Sample1"
    assert element_file.read == "R1"
    assert element_file.lane == ""  # Element files don't have lane info
    assert element_file.sample == ""  # Element files don't have sample info
    assert element_file.number == ""  # Element files don't have number info


# Tests for MergeGroup
def test_merge_group_init(illumina_pair):
    """Test initialization of MergeGroup."""
    r1_path, r2_path = illumina_pair
    r1_file = preprocessing.IlluminaFile(r1_path)
    r2_file = preprocessing.IlluminaFile(r2_path)

    merge_group = preprocessing.MergeGroup("Sample1", [r1_file, r2_file])
    assert merge_group.name == "Sample1"
    assert len(merge_group.files) == 2
    assert merge_group.merged_file is None


def test_merge_group_by_lane(illumina_pair):
    """Test _group_by_lane method of MergeGroup."""
    r1_path, r2_path = illumina_pair
    r1_file = preprocessing.IlluminaFile(r1_path)
    r2_file = preprocessing.IlluminaFile(r2_path)

    merge_group = preprocessing.MergeGroup("Sample1", [r1_file, r2_file])
    groups = merge_group._group_by_lane()

    # Should have one group for L001
    assert len(groups) == 1

    # Group should contain both R1 and R2
    assert len(groups[0]) == 2

    # Files should be in the group
    assert r1_file in groups[0]
    assert r2_file in groups[0]


def test_merge_group_merge(illumina_pair, merged_dir, log_dir):
    """Test merge method of MergeGroup."""
    r1_path, r2_path = illumina_pair
    r1_file = preprocessing.IlluminaFile(r1_path)
    r2_file = preprocessing.IlluminaFile(r2_path)

    merge_group = preprocessing.MergeGroup("Sample1", [r1_file, r2_file])
    result = merge_group.merge(
        merged_directory=merged_dir,
        log_directory=log_dir,
        format="fastq",
        algo="fastp",
        verbose=True,
    )

    # Check result
    expected_path = os.path.join(merged_dir, "Sample1.fastq")
    assert result == expected_path
    assert merge_group.merged_file == expected_path


# Tests for group_fastq_pairs
def test_group_fastq_pairs(multi_sample_fastqs):
    """Test grouping of FASTQ pairs."""
    r1_file = preprocessing.IlluminaFile(multi_sample_fastqs["r1_path"])
    r2_file = preprocessing.IlluminaFile(multi_sample_fastqs["r2_path"])
    r1_file2 = preprocessing.IlluminaFile(multi_sample_fastqs["r1_path2"])
    r2_file2 = preprocessing.IlluminaFile(multi_sample_fastqs["r2_path2"])

    files = [r1_file, r2_file, r1_file2, r2_file2]
    merge_groups = preprocessing.group_fastq_pairs(files)

    # Should have two merge groups
    assert len(merge_groups) == 2

    # Check that groups are correctly formed
    sample1_group = next((g for g in merge_groups if g.name == "Sample1"), None)
    sample2_group = next((g for g in merge_groups if g.name == "Sample2"), None)

    assert sample1_group is not None
    assert sample2_group is not None

    assert len(sample1_group.files) == 2
    assert len(sample2_group.files) == 2


def test_merge_fastqs(multi_sample_fastqs, merged_dir, log_dir):
    """Test merge_fastqs function."""

    result = preprocessing.merge_fastqs(
        files=multi_sample_fastqs["dirname"],
        output_directory=merged_dir,
        log_directory=log_dir,
        schema="illumina",
        verbose=True,
    )

    # Should return a list of merged files
    assert len(result) > 0
    # The merged file should be in the output directory
    assert os.path.dirname(result[0]) == merged_dir


def test_merge_fastqs_fastp(illumina_pair, merged_dir, log_dir):
    """Test merge_fastqs_fastp function."""
    r1_path, r2_path = illumina_pair

    merged_file = os.path.join(merged_dir, "merged.fastq")
    result = preprocessing.merge_fastqs_fastp(
        forward=r1_path,
        reverse=r2_path,
        merged=merged_file,
        log_directory=log_dir,
        name="Sample1",
    )

    # Check result
    assert result == merged_file


@mock.patch("abutils.tools.preprocessing.sp.Popen")
def test_merge_fastqs_vsearch(mock_popen, illumina_pair, merged_dir):
    """Test merge_fastqs_vsearch function."""
    r1_path, r2_path = illumina_pair

    # Need to mock subprocess.Popen for external commands
    mock_process = mock.MagicMock()
    mock_process.communicate.return_value = (b"", b"")
    mock_process.returncode = 0
    mock_popen.return_value = mock_process

    merged_file = os.path.join(merged_dir, "merged.fasta")
    result = preprocessing.merge_fastqs_vsearch(
        forward=r1_path,
        reverse=r2_path,
        merged_file=merged_file,
        output_format="fasta",
    )

    # Check result
    assert result == merged_file

    # Check that sp.Popen was called with vsearch command
    cmd = mock_popen.call_args[0][0]
    assert "vsearch" in cmd
    assert "--fastq_mergepairs" in cmd
    assert r1_path in cmd
    assert r2_path in cmd
    assert merged_file in cmd


# Test helper functions
def test_deduplicate_sequences():
    """Test _deduplicate_sequences function."""
    df = pl.DataFrame(
        {"sequence_id": ["seq1", "seq2", "seq3"], "sequence": ["ACGT", "ACGT", "TGCA"]}
    )

    # Test without keeping read numbers
    result = preprocessing._deduplicate_sequences(df, ["sequence"], False)
    assert result.height == 2  # Should have 2 unique sequences

    # Test with keeping read numbers
    result = preprocessing._deduplicate_sequences(df, ["sequence"], True)
    assert result.height == 2  # Should have 2 unique sequences
    assert "count" in result.columns  # Should have a count column
    assert result.select("count").to_series().to_list() == [2, 1]  # Check count values


# # Test for deduplicate and reduction functions
# @mock.patch("abutils.tools.preprocessing.to_fasta")
# def test_deduplicate(mock_to_fasta, airr_files, temp_dir):
#     """Test deduplicate function."""
#     # Mock to_fasta to avoid actual file writing
#     mock_to_fasta.return_value = None

#     output_dir = os.path.join(temp_dir, "dedup")

#     preprocessing.deduplicate(
#         # project_folder=airr_files["sample1_path"],  # Test file path input
#         project_folder=airr_files["dirname"],  # Test file path input
#         output="dedup",
#         output_format="fasta",
#         debug=True,
#     )

#     # Check that to_fasta was called with deduped data
#     assert mock_to_fasta.called
#     # The first argument should be a list of tuples with sequence_id and sequence
#     # We can't assert exact values because of potential mock issues, but we can check structure
#     args = mock_to_fasta.call_args[0]
#     assert isinstance(args[0], list)


@mock.patch("abutils.tools.cluster.cluster")
@mock.patch("abutils.tools.alignment.make_consensus")
def test_process_chains(mock_make_consensus, mock_cluster):
    """Test _process_chains function."""
    # Setup mock cluster
    mock_cluster_obj = mock.MagicMock()
    mock_cluster_obj.size = 3
    mock_cluster_obj.centroid = mock.MagicMock()
    mock_cluster_obj.centroid.id = "seq1"
    mock_cluster_obj.centroid.sequence = "ACGT"
    mock_cluster_obj.sequences = [
        Sequence("ACGT", id="seq1"),
        Sequence("ACGT", id="seq2"),
        Sequence("ACGT", id="seq3"),
    ]

    mock_cluster.return_value = [mock_cluster_obj]

    mock_consensus = Sequence("ACGT", id="consensus")
    mock_make_consensus.return_value = mock_consensus

    # Create test dataframe
    df = pl.DataFrame(
        {
            "sequence_id": ["seq1", "seq2", "seq3"],
            "sequence": ["ACGT", "ACGT", "ACGT"],
            "v_gene": ["IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01"],
            "j_gene": ["IGHJ1*01", "IGHJ1*01", "IGHJ1*01"],
            "vj_bin": [
                "IGHV1-1*01_IGHJ1*01",
                "IGHV1-1*01_IGHJ1*01",
                "IGHV1-1*01_IGHJ1*01",
            ],
        }
    )

    # Test centroid mode
    result = preprocessing._process_chains(
        "test",
        df,
        min_cluster_size=2,
        clustering_threshold=0.9,
        consentroid="centroid",
        cluster_sizes_separator="|",
        keep_cluster_sizes=True,
        output_format="fasta",
        debug=True,
    )

    assert len(result) == 1
    assert result[0].id == "seq1|3"
    assert result[0].sequence == "ACGT"

    # Test consensus mode
    result = preprocessing._process_chains(
        "test",
        df,
        min_cluster_size=2,
        clustering_threshold=0.9,
        consentroid="consensus",
        cluster_sizes_separator="|",
        keep_cluster_sizes=True,
        output_format="fasta",
        debug=True,
    )

    assert len(result) == 1
    assert result[0].id == "consensus|3"
    assert result[0].sequence == "ACGT"


def test_process_chain_group():
    """Test _process_chain_group function."""
    # Create test dataframe
    df = pl.DataFrame(
        {
            "sequence_id": ["seq1", "seq2", "seq3"],
            "sequence": ["ACGT", "ACGT", "TGCA"],
            "v_gene": ["IGHV1-1*01", "IGHV1-1*01", "IGHV2-1*01"],
            "j_gene": ["IGHJ1*01", "IGHJ1*01", "IGHJ2*01"],
            "locus": ["IGH", "IGH", "IGH"],
        }
    )

    with mock.patch(
        "abutils.tools.preprocessing._process_chains"
    ) as mock_process_chains:
        mock_process_chains.return_value = [Sequence("ACGT", id="seq1")]

        # Test without UMI
        result = preprocessing._process_chain_group(
            df,
            chain_type="Heavy",
            umi=False,
            min_cluster_size=2,
            clustering_threshold=0.9,
            consentroid="centroid",
            cluster_sizes_separator="|",
            keep_cluster_sizes=False,
            output_format="fasta",
            debug=True,
        )

        assert result == [Sequence("ACGT", id="seq1")]

        # Check that vj_bin was created correctly
        assert "vj_bin" in mock_process_chains.call_args[0][1].columns
        vj_bins = (
            mock_process_chains.call_args[0][1].select("vj_bin").to_series().to_list()
        )
        assert "IGHV1-1*01_IGHJ1*01" in vj_bins
        assert "IGHV2-1*01_IGHJ2*01" in vj_bins

        # Test with UMI
        df_with_umi = df.with_columns(pl.lit("UMI1").alias("umi"))

        result = preprocessing._process_chain_group(
            df_with_umi,
            chain_type="Heavy",
            umi=True,
            min_cluster_size=2,
            clustering_threshold=0.9,
            consentroid="centroid",
            cluster_sizes_separator="|",
            keep_cluster_sizes=False,
            output_format="fasta",
            debug=True,
        )

        assert result == [Sequence("ACGT", id="seq1")]

        # Check that vj_bin includes UMI
        vj_bins = (
            mock_process_chains.call_args[0][1].select("vj_bin").to_series().to_list()
        )
        assert any("UMI1" in bin_name for bin_name in vj_bins)


# @mock.patch("abutils.tools.preprocessing.to_fasta")
# @mock.patch("abutils.tools.preprocessing._process_chain_group")
# def test_reduction(mock_process_chain_group, mock_to_fasta, airr_files, temp_dir):
#     """Test reduction function."""
#     mock_process_chain_group.return_value = [
#         Sequence("ACGT", id="seq1"),
#         Sequence("TGCA", id="seq2"),
#     ]

#     mock_to_fasta.return_value = None

#     output_dir = os.path.join(temp_dir, "reduced")

#     preprocessing.reduction(
#         project_folder=airr_files["sample1_path"],  # Test file path input
#         output="reduced",
#         output_format="fasta",
#         debug=True,
#     )

#     # Check that _process_chain_group and to_fasta were called
#     assert mock_process_chain_group.called
#     assert mock_to_fasta.called
