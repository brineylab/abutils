#!/usr/bin/python
# filename: preprocessing.py

#
# Copyright (c) 2024 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


import datetime
import os
import subprocess as sp
import tempfile
import time
from typing import Iterable, Optional, Union

import polars as pl
from natsort import natsorted
from tqdm.auto import tqdm

from .. import Sequence
from ..bin import get_path as get_binary_path
from ..io import (
    concatenate_files,
    delete_files,
    from_polars,
    list_files,
    make_dir,
    rename_file,
    to_fasta,
    to_polars,
)
from ..tools import alignment, cluster

__all__ = ["deduplicate", "reduction", "merge_fastqs", "group_fastq_pairs"]


class FASTQFile:
    """
    Base class for processing paired FASTQ filenames.

    Parameters
    ----------
    file : str
        Path to a FASTQ file.

    """

    def __init__(self, file: str):
        self.file = file
        self.path = os.path.abspath(file)
        self.basename = os.path.basename(self.path)
        self.dir = os.path.dirname(self.path)
        self.filename = self.basename.rstrip(".gz").rstrip(".fastq").rstrip(".fq")

    def __eq__(self, other):
        return self.name == other.name


class IlluminaFile(FASTQFile):
    """
    Class for processing Illumina-formatted FASTQ filenames.
    Also compatible with Element's "legacy" filename structure.

    Parameters
    ----------
    file : str
        Path to a FASTQ file.

    """

    def __init__(self, file: str):
        super().__init__(file)
        self.schema = "illumina"

    @property
    def name(self):
        return "_".join(self.filename.split("_")[:-4])

    @property
    def number(self):
        return self.filename.split(".")[0].split("_")[-1]

    @property
    def read(self):
        return self.filename.split("_")[-2]

    @property
    def lane(self):
        return self.filename.split("_")[-3]

    @property
    def sample(self):
        return self.filename.split("_")[-4]


class ElementFile(FASTQFile):
    """
    Class for processing Element-formatted FASTQ filenames.

    Parameters
    ----------
    file : str
        Path to a FASTQ file.

    """

    def __init__(self, file: str):
        super().__init__(file)
        self.schema = "element"

    @property
    def name(self):
        return "_".join(self.filename.split("_")[:-1])

    @property
    def number(self):
        return ""

    @property
    def read(self):
        return self.filename.split("_")[-1]

    @property
    def lane(self):
        return ""

    @property
    def sample(self):
        return ""


class MergeGroup:
    """
    Group of FASTQ files to be merged.

    Parameters
    ----------
    name : str
        Name of the sample.

    files : Iterable[FASTQFile]
        Iterable of FASTQFile objects.

    """

    def __init__(self, name: str, files: Iterable[FASTQFile]):
        self.name = name
        self.files = files
        self.merged_file = None  # assigned by merge function

    def merge(
        self,
        merged_directory: str,
        log_directory: Optional[str] = None,
        format: str = "fastq",
        algo: str = "fastp",
        binary_path: Optional[str] = None,
        minimum_overlap: int = 30,
        allowed_mismatches: int = 5,
        allowed_mismatch_percent: float = 20.0,
        trim_adapters: bool = True,
        adapter_file: Optional[str] = None,
        quality_trim: bool = True,
        window_size: int = 4,
        quality_cutoff: int = 20,
        merge_args: Optional[str] = None,
        compress_output: bool = False,
        verbose: bool = False,
        debug: bool = False,
    ) -> str:
        groups = self._group_by_lane()
        n_groups = len(groups)
        if verbose:
            print(f"\n{self.name}")
            # print("-" * len(self.name))
            if n_groups > 1:
                groups = tqdm(
                    groups,
                    desc="  - merging lanes",
                    total=n_groups,
                    unit="lane",
                    leave=False,
                )
            else:
                print("  - merging paired FASTQs")

        # merge function
        if algo.lower() == "fastp":
            merge_func = merge_fastqs_fastp
        else:
            raise ValueError(f"Invalid merge algorithm: {algo}. Must be 'fastp'.")

        # merge files
        merged_files = []
        compress_suffix = ".gz" if compress_output else ""
        self.merged_file = os.path.join(
            merged_directory, f"{self.name}.{format.lower()}{compress_suffix}"
        )
        for group in groups:
            r1, r2 = natsorted(group, key=lambda x: x.read)
            # add name and compression suffix, if necessary
            name = self.name if n_groups == 1 else f"{self.name}_{r1.lane}"
            merged_path = os.path.join(
                merged_directory, f"{name}.{format.lower()}{compress_suffix}"
            )
            # merge
            merge_func(
                r1.path,
                r2.path,
                merged_path,
                name=name,
                log_directory=log_directory,
                trim_adapters=trim_adapters,
                quality_trim=quality_trim,
                window_size=window_size,
                quality_cutoff=quality_cutoff,
                minimum_overlap=minimum_overlap,
                allowed_mismatches=allowed_mismatches,
                allowed_mismatch_percent=allowed_mismatch_percent,
                adapter_file=adapter_file,
                binary_path=binary_path,
                additional_args=merge_args,
                debug=debug,
            )
            merged_files.append(merged_path)

        # concatenate merged files if dataset has multiple lanes
        if len(merged_files) > 1:
            if verbose:
                print("  - concatenating merged files")
            concatenate_files(merged_files, self.merged_file)
            delete_files(merged_files)
        else:
            rename_file(merged_files[0], self.merged_file)
        return self.merged_file

    def _group_by_lane(self):
        lane_dict = {}
        for f in self.files:
            if f.lane not in lane_dict:
                lane_dict[f.lane] = []
            lane_dict[f.lane].append(f)
        return [
            group for lane, group in natsorted(lane_dict.items(), key=lambda x: x[0])
        ]


def merge_fastqs(
    files: Union[str, Iterable],
    output_directory: str,
    output_format: str = "fastq",
    log_directory: Optional[str] = None,
    schema: str = "illumina",
    algo: str = "fastp",
    binary_path: Optional[str] = None,
    merge_args: Optional[str] = None,
    minimum_overlap: int = 30,
    allowed_mismatches: int = 5,
    allowed_mismatch_percent: float = 20.0,
    trim_adapters: bool = True,
    adapter_file: Optional[str] = None,
    quality_trim: bool = True,
    window_size: int = 4,
    quality_cutoff: int = 20,
    compress_output: bool = False,
    debug: bool = False,
    verbose: bool = False,
) -> Iterable[str]:
    """
    Merge paired-end fastq files.

    Parameters
    ----------
    files : Union[str, Iterable]
        Path to a directory containing paired-end fastq files, or an iterable of paired-end fastq files.

    output_directory : str
        Path to the directory in which to save the merged files.

    output_format : str, optional
        Output format. Must be 'fastq' or 'fasta'. Default is 'fastq'.

    log_directory : str, optional
        Path to the directory in which to save the fastp reports.
        If not provided, fastp reports will not be saved.

    schema : str, optional
        Schema of the file names. Can be 'illumina' or 'element'. Default is 'illumina'.

    algo : str, optional
        Algorithm to use for merging. Must be 'fastp'. Default is 'fastp'.

    binary_path : str, optional
        Path to the merge algorithm binary. If not provided, the binary packaged with abutils will be used.

    merge_args : str, optional
        Additional arguments (as a string) to pass to the merge function.

    minimum_overlap : int, optional
        Minimum overlap between reads. Default is 30.

    allowed_mismatches : int, optional
        Allowed mismatches between reads. Default is 5.

    allowed_mismatch_percent : float, optional
        Allowed mismatch percentage between reads. Default is 20.0.

    trim_adapters : bool, optional
        If True, trim adapters. Default is True.

    adapter_file : str, optional
        Path to a FASTA file containing adapter sequences. If not provided, the default
        adapter sequences (Illumina TruSeq) will be used.

    quality_trim : bool, optional
        If True, trim low-quality bases. Default is True.

    window_size : int, optional
        Sliding window size for quality trimming. Default is 4.

    quality_cutoff : int, optional
        Mean quality cutoff for quality trimming. Default is 20.

    compress_output : bool, optional
        If True, compress the output file using gzip. Default is False.

    debug : bool, optional
        If True, print debug output. Default is False.

    verbose : bool, optional
        If True, print verbose output. Default is False.

    Returns
    -------
    list
        A list of paths to the merged FASTA/Q files.

    """
    # process input files
    if isinstance(files, str):
        if not os.path.isdir(files):
            err = f"The supplied file path ({files}) does not exist or is not a directory."
            raise ValueError(err)
        files = list_files(files, extension=[".fastq", ".fq", ".fastq.gz", ".fq.gz"])
    if schema.lower() == "illumina":
        files = [IlluminaFile(f) for f in files]
    elif schema.lower() == "element":
        files = [ElementFile(f) for f in files]
    else:
        raise ValueError(f"Invalid schema: {schema}. Must be 'illumina' or 'element'.")
    # group files by sample
    file_pairs = group_fastq_pairs(files, verbose=verbose)
    if verbose:
        print("--------------------")
        print("    MERGE FASTQs    ")
        print("--------------------")
        file_pairs = tqdm(file_pairs, leave=True)
    # merge files
    make_dir(output_directory)
    merged_files = []
    for fp in file_pairs:
        merged_file = fp.merge(
            merged_directory=output_directory,
            format=output_format,
            algo=algo,
            binary_path=binary_path,
            log_directory=log_directory,
            minimum_overlap=minimum_overlap,
            allowed_mismatches=allowed_mismatches,
            allowed_mismatch_percent=allowed_mismatch_percent,
            trim_adapters=trim_adapters,
            adapter_file=adapter_file,
            quality_trim=quality_trim,
            window_size=window_size,
            quality_cutoff=quality_cutoff,
            merge_args=merge_args,
            compress_output=compress_output,
            debug=debug,
            verbose=verbose,
        )
        merged_files.append(merged_file)
    return merged_files


def group_fastq_pairs(
    files: Iterable[FASTQFile], verbose: bool = False
) -> Iterable[MergeGroup]:
    """
    Group paired-end fastq files by sample name. If a sample has multiple lanes, the files will be combined.

    Parameters
    ----------
    files : Union[str, Iterable]
        Path to a directory containing paired-end fastq files, or an iterable of paired-end fastq files.

    verbose : bool, optional
        If True, print verbose output. Default is False.

    Returns
    -------
    list
        A list of MergeGroup objects.

    """
    group_dict = {}
    for f in files:
        if f.name not in group_dict:
            group_dict[f.name] = []
        group_dict[f.name].append(f)
    groups = [MergeGroup(name, group_files) for name, group_files in group_dict.items()]
    return natsorted(groups, key=lambda x: x.name)


def merge_fastqs_vsearch(
    forward: str,
    reverse: str,
    merged_file: str,
    output_format: str = "fasta",
    binary_path: Optional[str] = None,
    minimum_overlap: int = 30,
    allowed_mismatches: int = 5,
    allowed_mismatch_percent: float = 20.0,
    additional_args: Optional[str] = None,
    debug: bool = False,
) -> str:
    """
    Merge paired-end reads using vsearch.

    Parameters
    ----------
    forward : str
        Path to the forward read file.

    reverse : str
        Path to the reverse read file.

    merged_file : str
        Path to the merged read file.

    output_format : str, optional
        Output format. Must be 'fasta' or 'fastq'. Default is 'fasta'.

    binary_path : str, optional
        Path to a vsearch binary. If not provided, the vsearch binary packaged with abutils will be used.

    additional_args : str, optional
        Additional arguments (as a string) to pass to vsearch.

    debug : bool, optional
        If True, print vsearch stdout and stderr. Default is False.

    Returns
    -------
    str
        Path to the merged read file.

    """
    # validate input files
    if not os.path.isfile(forward):
        err = f"The supplied forward read file path ({forward}) does not exist or is not a file."
        raise ValueError(err)
    if not os.path.isfile(reverse):
        err = f"The supplied reverse read file path ({reverse}) does not exist or is not a file."
        raise ValueError(err)
    # make output directory
    out_dir = os.path.dirname(merged_file)
    make_dir(out_dir)
    # get the vsearch binary
    if binary_path is None:
        binary_path = get_binary_path("vsearch")
    # compile the vsearch command
    cmd = f"{binary_path} --fastq_mergepairs {forward} --reverse {reverse}"
    if output_format.lower() == "fasta":
        cmd += f" --fastaout {merged_file}"
    elif output_format.lower() == "fastq":
        cmd += f" --fastqout {merged_file}"
    else:
        err = f"Invalid output format: {output_format}. Must be 'fasta' or 'fastq'."
        raise ValueError(err)
    if additional_args is not None:
        cmd += f" {additional_args}"
    # merge reads
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(stdout.decode())
        print(stderr.decode())
    if p.returncode != 0:
        err = f"Error merging reads with vsearch: {stderr.decode()}"
        raise ValueError(err)
    return merged_file


def merge_fastqs_fastp(
    forward: str,
    reverse: str,
    merged: str,
    binary_path: Optional[str] = None,
    minimum_overlap: int = 30,
    allowed_mismatches: int = 5,
    allowed_mismatch_percent: int = 20,
    correct_overlap_region: bool = False,
    trim_adapters: bool = True,
    adapter_file: Optional[str] = None,
    quality_trim: bool = True,
    window_size: int = 4,
    quality_cutoff: int = 20,
    name: Optional[str] = None,
    log_directory: Optional[str] = None,
    additional_args: Optional[str] = None,
    debug: bool = False,
) -> str:
    """
    Merge paired-end reads using vsearch.

    Parameters
    ----------
    forward : str
        Path to the forward read file.

    reverse : str
        Path to the reverse read file.

    merged : str
        Path to the merged read file. If the merged file ends in ".gz", the merged file
        will be gzip compressed.

    binary_path : str, optional
        Path to a fastp binary. If not provided, the fastp binary packaged with abutils will be used.

    minimum_overlap : int, optional
        Minimum overlap between reads. Default is 30.

    allowed_mismatches : int, optional
        Allowed mismatches between reads. Default is 5.

    allowed_mismatch_percent : float, optional
        Allowed mismatch percentage between reads. Default is 20.0.

    correct_overlap_region : bool, optional
        If True, correct the overlap region. Default is False.

    trim_adapters : bool, optional
        If True, trim adapters. Default is True.

    adapter_file : str, optional
        Path to a FASTA file containing adapter sequences. If not provided, the default
        adapter sequences (Illumina TruSeq) will be used.

    quality_trim : bool, optional
        If True, trim low-quality bases. Default is True.

    window_size : int, optional
        Sliding window size for quality trimming. Default is 4.

    quality_cutoff : int, optional
        Mean quality cutoff for quality trimming. Default is 20.

    name : str, optional
        Name of the sample. If not provided, the name of the merged file will be used.

    log_directory : str, optional
        Path to the directory in which to save the fastp report.
        If not provided, fastp reports will not be saved.

    additional_args : str, optional
        Additional arguments (as a string) to pass to fastp.

    debug : bool, optional
        If True, print vsearch stdout and stderr. Default is False.

    Returns
    -------
    str
        Path to the merged read file.

    """
    # validate input files
    if not os.path.isfile(forward):
        err = f"The supplied forward read file path ({forward}) does not exist or is not a file."
        raise ValueError(err)
    if not os.path.isfile(reverse):
        err = f"The supplied reverse read file path ({reverse}) does not exist or is not a file."
        raise ValueError(err)

    # make output directory
    merged = os.path.abspath(merged)
    out_dir = os.path.dirname(merged)
    make_dir(out_dir)

    # get the vsearch binary
    if binary_path is None:
        binary_path = get_binary_path("fastp")

    # compile the vsearch command
    cmd = f"{binary_path} -i '{forward}' -I '{reverse}' --merge --merged_out '{merged}'"
    cmd += f" --overlap_len_require {int(minimum_overlap)}"
    cmd += f" --overlap_diff_limit {int(allowed_mismatches)}"
    cmd += f" --overlap_diff_percent_limit {int(allowed_mismatch_percent)}"
    if correct_overlap_region:
        cmd += " --correction"

    # adapters
    if trim_adapters:
        if adapter_file is None:
            cmd += " --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
            cmd += " --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        else:
            cmd += f" --adapter_fasta {adapter_file}"
    else:
        cmd += " --disable_adapter_trimming"

    # quality
    if quality_trim:
        cmd += " --cut_tail"
        cmd += f" --cut_tail_window_size {int(window_size)}"
        cmd += f" --cut_tail_mean_quality {int(quality_cutoff)}"
    else:
        cmd += " --disable_quality_filtering"

    # log
    if log_directory is None:
        log_directory = tempfile.mkdtemp()
    make_dir(log_directory)
    if name is None:
        if merged.endswith(".gz"):
            name = ".".join(os.path.basename(merged).split(".")[:-2])
        else:
            name = ".".join(os.path.basename(merged).split(".")[:-1])
    cmd += f" --html '{os.path.join(log_directory, name)}_fastp-report.html'"
    cmd += f" --json '{os.path.join(log_directory, name)}_fastp-report.json'"

    # additional CLI args
    if additional_args is not None:
        cmd += f" {additional_args}"

    # merge reads
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if debug:
        print(stdout.decode())
        print(stderr.decode())
    if p.returncode != 0:
        err = f"Error merging reads with fastp: {stderr.decode()}"
        raise ValueError(err)
    return merged


#########################################
#                                       #
#         Helper functions              #
#                                       #
#########################################


def _deduplicate_sequences(
    df: pl.DataFrame, group_by_cols: list, keep_read_numbers: bool
):
    if not keep_read_numbers:
        return df.unique(subset=group_by_cols).sort("sequence")
    return (
        df.group_by(group_by_cols)
        .agg(
            pl.count().alias("count"),
            pl.col("sequence_id").first(),
        )
        .sort("sequence")
    )


def _write_fasta(
    df: pl.DataFrame,
    file_path: str,
    keep_read_numbers: bool,
    read_number_separator: str,
):
    if keep_read_numbers:
        df = df.with_columns(
            (
                pl.col("sequence_id")
                + read_number_separator
                + pl.col("count").cast(pl.Utf8)
            ).alias("sequence_id")
        )
    tuples_list = list(zip(df["sequence_id"].to_list(), df["sequence"].to_list()))
    to_fasta(tuples_list, file_path)
    print(f"Output written to FASTA file: {file_path}\n")


def _process_chains(
    step: str,
    chains: pl.DataFrame,
    min_cluster_size: int,
    clustering_threshold: float,
    consentroid: str,
    cluster_sizes_separator: str,
    keep_cluster_sizes: bool,
    output_format: str,
    debug: bool,
):
    """Internal function to generate a consentroid from a bin of sequences"""

    consentroids = []
    for vj in tqdm(
        chains["vj_bin"].unique().to_list(), desc=f"{step} chains processing"
    ):
        group = chains.filter(pl.col("vj_bin") == vj)
        if group.height < min_cluster_size:
            continue

        group_sequences = from_polars(group)
        clusters = cluster.cluster(group_sequences, threshold=clustering_threshold)
        for c in clusters:
            if c.size < min_cluster_size:
                continue
            if consentroid == "centroid":
                seq_id = (
                    f"{c.centroid.id}{cluster_sizes_separator}{c.size}"
                    if keep_cluster_sizes
                    else c.centroid.id
                )
                consentroids.append(Sequence(c.centroid.sequence, id=seq_id))
            elif consentroid == "consensus":
                _consensus = alignment.make_consensus(
                    [s for s in c.sequences],
                    algo="mafft",
                    alignment_kwargs={"mafft_bin": "mafft", "threads": 1},
                    debug=debug,
                )
                seq_id = (
                    f"{_consensus.id}{cluster_sizes_separator}{c.size}"
                    if keep_cluster_sizes
                    else _consensus.id
                )
                consentroids.append(Sequence(_consensus.sequence, id=seq_id))

        if debug:
            print(
                f"\tProcessed {len(clusters)} clusters for {vj}, corresponding to {len(group_sequences)} sequences"
            )

    if debug:
        print(f"\tTotal concentroids generated: {len(consentroids)}\n")

    if output_format == "fasta":
        return consentroids

    elif output_format == "airr":
        df = to_polars(consentroids)
        if keep_cluster_sizes and len(df) != 0:
            df = (
                df.with_columns(
                    pl.col("sequence_id")
                    .str.split_exact(
                        cluster_sizes_separator, 1
                    )  # Split into a struct with two fields
                    .alias("split_id")
                )
                .with_columns(
                    pl.col("split_id")
                    .struct.field("field_0")
                    .alias("sequence_id"),  # Extract the first field
                    pl.col("split_id")
                    .struct.field("field_1")
                    .cast(pl.Utf8)
                    .alias("count"),  # Extract the second field
                )
                .drop("split_id")
            )  # Drop the intermediate "split_id" column
        return df


def _process_chain_group(
    chains: pl.DataFrame,
    chain_type: str,
    umi: bool,
    min_cluster_size: int,
    clustering_threshold: float,
    consentroid: str,
    cluster_sizes_separator: str,
    keep_cluster_sizes: bool,
    output_format: str,
    debug: bool,
):
    # Generate vj_bin column
    if umi:
        chains = chains.with_columns(
            (pl.col("v_gene") + "_" + pl.col("j_gene") + "_" + pl.col("umi")).alias(
                "vj_bin"
            )
        )
    else:
        chains = chains.with_columns(
            (pl.col("v_gene") + "_" + pl.col("j_gene")).alias("vj_bin")
        )

    return _process_chains(
        chain_type,
        chains,
        min_cluster_size,
        clustering_threshold,
        consentroid,
        cluster_sizes_separator,
        keep_cluster_sizes,
        output_format,
        debug,
    )


#########################################
#                                       #
#           Main functions              #
#                                       #
#########################################


def deduplicate(
    project_folder: str,
    output: Optional[str] = None,
    output_format: str = "fasta",
    pool: bool = True,
    umi: bool = False,
    keep_read_numbers: bool = False,
    read_number_separator: str = "|",
    large_files: bool = False,
    debug: bool = False,
) -> None:
    """
    A polyvalent tool for deduplication of assigned reads.
    This function takes as input the AIRR-compliant
    tables and is specifically designed to handle extremely
    large files with a minimal footprint.

    Parameters
    ----------
    project_folder : str
        Path to the project folder containing AIRR-compliant tables.

    output : str, optional
        Subdirectory for output files. Created if non-existent. Defaults to None.

    output_format : str
        Either "fasta" or "airr". Default is "fasta".

    pool : bool, optional
        If True, pool all samples together. Defaults to True.

    umi : bool, optional
        If True, use UMI for deduplication. Defaults to False.

    keep_read_numbers : bool, optional
        If True, read numbers will be added to sequence names. Defaults to False.

    read_number_separator : str, optional
        Separator for read numbers in sequence names. Defaults to "|".

    large_files : bool, optional
        If True, optimize for large files (>100Go). Defaults to False.

    debug : bool, optional
        If True, print debug information. Defaults to False.

    Returns
    -------
    None. A FASTA- or TSV-formatted file is written to disk.
    """

    start = time.time()

    # Assert that output_format is valid
    if output_format not in ["fasta", "airr"]:
        raise ValueError("Invalid output_format. Must be 'fasta' or 'airr'.")

    # Preparing input and output file(s) / folder(s)
    if os.path.isfile(project_folder):
        project_folder = os.path.dirname(project_folder)
        if debug:
            print(f"Project folder is set to: {project_folder}")

    files = sorted(
        [
            f
            for f in list_files(
                project_folder,
                extension="tsv",
                recursive=True,
            )
            if "tmp" not in f
        ]
    )

    if len(files) == 0:
        print("No files found. Exiting.")
    elif len(files) == 1:
        pool = False

    if output:
        project_folder = os.path.join(project_folder, output)
        make_dir(project_folder)

    # Running main deduplication process
    total_sequences = 0
    pooled = []

    for file in files:
        sample = os.path.basename(file).split(".")[0]

        print(f"Processing {sample}")
        print("-" * (len(sample) + 11) + "\n")

        # Loading files
        if umi:
            try:
                df = pl.read_csv(
                    file,
                    columns=["sequence_id", "sequence", "umi"],
                    separator="\t",
                    null_values="None",
                    low_memory=True if large_files else False,
                )
            except Exception as e:
                print(f"Error reading file {file}: {e}")
                continue
        else:
            try:
                df = pl.read_csv(
                    file,
                    columns=["sequence_id", "sequence"],
                    separator="\t",
                    null_values="None",
                    low_memory=True if large_files else False,
                )
            except Exception as e:
                print(f"Error reading file {file}: {e}")
                continue

        print(f"Loaded {df.height:,} annotated sequences")
        total_sequences += df.shape[0]

        # Processing with deduplication in the presence of UMIs
        if umi:
            if debug:
                print(
                    f"Deduplicating with UMI ({df.unique(subset=['umi']).height} unique UMIs)"
                )

            df_unique = _deduplicate_sequences(
                df, ["sequence", "umi"], keep_read_numbers
            )

        # Processing with deduplication in the absence of UMIs
        else:
            df_unique = _deduplicate_sequences(df, ["sequence"], keep_read_numbers)

        msg = f"Found {df_unique.height:,} unique sequences"
        if pool:
            msg += " added to the pool\n"
        print(msg)

        if pool:
            pooled.append(df_unique)
        else:
            # Saving deduplicated single output to FASTA file
            if output_format == "fasta":
                fasta_file = os.path.join(project_folder, sample) + ".fasta"

                _write_fasta(
                    df_unique, fasta_file, keep_read_numbers, read_number_separator
                )

            # Saving deduplicated single output to TSV file
            elif output_format == "airr":
                tsv_file = os.path.join(project_folder, sample) + ".tsv"
                # Convert df_unique to LazyFrame for more efficient processing
                df_unique_lazy = pl.LazyFrame(df_unique)
                # Scan the complete original AIRR table
                new_df = pl.scan_csv(
                    file, separator="\t", low_memory=True if large_files else False
                )

                # Filter new_df to only include rows with matching sequence_ids
                filtered_df = new_df.join(
                    df_unique_lazy.select(["sequence_id"]), on="sequence_id", how="semi"
                ).select(new_df.columns)

                # Add the 'count' column to the filtered DataFrame (only for matching rows)
                if keep_read_numbers:
                    filtered_df = (
                        filtered_df.join(
                            df_unique_lazy.select(["sequence_id", "count"]),
                            on="sequence_id",
                            how="left",
                        )
                        .with_columns(
                            pl.col("count").fill_null(
                                0
                            )  # Replace null values with 0 if necessary
                        )
                        .rename(
                            {"count": "duplicates"}  # Rename the column
                        )
                    )

                # Sink the resulting DataFrame to a TSV file
                filtered_df.collect().write_csv(tsv_file, separator="\t")
                print(f"Output written to TSV file: {tsv_file}\n")

    if pool:
        pool_df = pl.concat(pooled)

        if not keep_read_numbers:
            pool_unique = pool_df.unique(subset=["sequence"]).sort("sequence")
        else:
            pool_unique = (
                pool_df.group_by("sequence")
                .agg(
                    pl.count().alias("count"),
                    pl.col("sequence_id").first(),
                )
                .sort("sequence")
            )

        print(f"\nFound {pool_unique.height:,} unique sequences in pooled data")

        # Saving deduplicated pooled output to FASTA file
        if output_format == "fasta":
            fasta_file = os.path.join(project_folder, "deduplicated_pool.fasta")
            _write_fasta(
                pool_unique, fasta_file, keep_read_numbers, read_number_separator
            )

        # Saving deduplicated pooled output to TSV file
        elif output_format == "airr":
            tsv_file = os.path.join(project_folder, "deduplicated_pool.tsv")

            # Convert df_unique to LazyFrame for more efficient processing
            df_unique_lazy = pl.LazyFrame(pool_unique)

            # Scan and concatenate multiple AIRR tables
            new_df = pl.concat(
                [
                    pl.scan_csv(
                        file, separator="\t", low_memory=True if large_files else False
                    )
                    for file in files
                ]
            )

            # Filter new_df to only include rows with matching sequence_ids
            filtered_df = new_df.join(
                df_unique_lazy.select(["sequence_id"]), on="sequence_id", how="semi"
            ).select(new_df.columns)

            # Add the 'count' column to the filtered DataFrame (only for matching rows)
            if keep_read_numbers:
                filtered_df = (
                    filtered_df.join(
                        df_unique_lazy.select(["sequence_id", "count"]),
                        on="sequence_id",
                        how="left",
                    )
                    .with_columns(
                        pl.col("count").fill_null(
                            0
                        )  # Replace null values with 0 if necessary
                    )
                    .rename(
                        {"count": "duplicates"}  # Rename the column
                    )
                )

            # Sink the resulting DataFrame to a TSV file
            filtered_df.collect().write_csv(tsv_file, separator="\t")
            print(f"Output written to TSV file: {tsv_file}\n")

    # Finalizing
    end = time.time()
    duration = end - start
    print(
        f"Deduplication complete. Total elapsed time: {datetime.timedelta(seconds=int(duration))}."
    )

    return


def reduction(
    project_folder: str,
    output: Optional[str] = None,
    output_format: str = "fasta",
    pool: bool = True,
    umi: bool = False,
    keep_cluster_sizes: bool = False,
    cluster_sizes_separator: str = "|",
    min_cluster_size: int = 3,
    clustering_threshold: float = 0.975,
    consentroid: str = "centroid",
    large_files: bool = False,
    debug: bool = False,
) -> None:
    """
    This function takes as an input AIRR-compliant tables (tsv) and proceeds to
    data reduction by clustering sequences to a high identity threshold.

    This is specifically designed to handle large files with minimal footprint.
    Preclustering can be applied to increase performance over large datasets.

    Parameters:
    project_folder (str): Path to the project folder containing AIRR-compliant tables.
    output (str, optional): Subdirectory for output files. Created if non-existent. Defaults to None.
    output_format (str): Either "fasta" or "airr". Default is "fasta".
    pool (bool, optional): If True, pool all samples together. Defaults to True.
    umi (bool, optional): If True, use UMI for clustering. Defaults to False.
    keep_cluster_sizes (bool, optional): If True, cluster sizes will be added to sequence names. Defaults to False.
    cluster_sizes_separator (str, optional): Separator for cluster sizes in sequence names. Defaults to '|'.
    min_cluster_size (int, optional): Minimum cluster size to consider. Defaults to 3.
    clustering_threshold (float, optional): Identity threshold for clustering. Defaults to 0.975.
    consentroid (str, optional): Method to determine cluster representative ('centroid' or 'consensus'). Defaults to 'centroid'.
    large_files (bool, optional): If True, optimize for large files (>100Go). Defaults to False.
    debug (bool, optional): If True, print debug information. Defaults to False.

    Returns: None. A fasta file is written to disk.
    """

    start = time.time()

    if umi:
        cluster_sizes_separator += "umi_count="
    else:
        cluster_sizes_separator += "cluster_size="

    # Assert that output_format is valid
    if output_format not in ["fasta", "airr"]:
        raise ValueError("Invalid output_format. Must be 'fasta' or 'airr'.")

    # Assert that clustering_threshold is valid
    if not (0 < clustering_threshold <= 1):
        raise ValueError("Clustering_threshold must be between 0 and 1.")

    # Assert that consentroid is valid
    if consentroid not in ["centroid", "consensus"]:
        raise ValueError("Consentroid must be either 'centroid' or 'consensus'.")

    # Check for incompatible consentroid and output_format arguments
    if consentroid == "consensus" and output_format == "airr":
        raise ValueError(
            "The 'consensus' consentroid method is incompatible with the 'airr' output format. "
            "Please use 'centroid' or choose 'fasta' as the output format."
        )

    # Preparing input and output file(s) / folder(s)
    if os.path.isfile(project_folder):
        project_folder = os.path.dirname(project_folder)
        if debug:
            print(f"Project folder is set to: {project_folder}")

    files = sorted(
        [
            f
            for f in list_files(
                project_folder,
                extension="tsv",
                recursive=True,
            )
            if not "tmp" in f
        ]
    )

    if len(files) == 0:
        print("No files found. Exiting.")
    elif len(files) == 1:
        pool = False

    if output:
        project_folder = os.path.join(project_folder, output)
        make_dir(project_folder)

    total_sequences = 0

    if pool:
        pooled_heavies = []
        pooled_lights = []

    # Processing files individually...
    for file in files:
        sample = os.path.basename(file).split(".")[0]

        print(f"Processing {sample}")
        print("-" * (len(sample) + 11) + "\n")

        keys = ["sequence_id", "v_gene", "j_gene", "locus", "sequence"]
        if umi:
            keys.append("umi")

        df = pl.read_csv(
            file,
            columns=keys,
            separator="\t",
            null_values="None",
            low_memory=True if large_files else False,
        )

        if not all(col in df.columns for col in keys):
            print(f"File {file} is missing required columns. Skipping...")
            continue

        print(f"Loaded {df.height:,} annotated sequences")
        total_sequences += df.height

        heavies = df.filter(pl.col("locus") == "IGH")
        lights = df.filter(pl.col("locus") != "IGH")

        if pool:
            pooled_heavies.append(heavies)
            pooled_lights.append(lights)

        else:
            sample_consentroids = []

            # Process individual heavy chains
            heavy_consentroids = _process_chain_group(
                heavies,
                "Heavy",
                umi,
                min_cluster_size,
                clustering_threshold,
                consentroid,
                cluster_sizes_separator,
                keep_cluster_sizes,
                output_format,
                debug,
            )
            sample_consentroids.extend(heavy_consentroids)

            # Process individual light chains
            light_consentroids = _process_chain_group(
                lights,
                "Lights",
                umi,
                min_cluster_size,
                clustering_threshold,
                consentroid,
                cluster_sizes_separator,
                keep_cluster_sizes,
                output_format,
                debug,
            )
            sample_consentroids.extend(light_consentroids)

            # Skipping file creation if no concentroids have been generated
            if not sample_consentroids:
                print(f"No clusters found for sample {sample}. Skipping output.\n")
                continue

            if output_format == "fasta":
                # Save consentroids to FASTA file
                fasta_file = os.path.join(project_folder, sample) + "_reduced.fasta"
                to_fasta(sample_consentroids, fasta_file)

            elif output_format == "airr":
                # Saving centroids to TSV file (consensus can't be exported as AIRR-annotated TSV table as they haven't been annotated yet)
                tsv_file = os.path.join(project_folder, sample + "_reduced.tsv")

                # Convert centroids to LazyFrame for more efficient processing
                centroids_lazy = pl.LazyFrame(
                    sample_consentroids,
                )

                # Scan the complete original AIRR table
                new_df = pl.scan_csv(
                    file, separator="\t", low_memory=True if large_files else False
                )

                # Filter new_df to only include rows with matching sequence_ids
                filtered_df = new_df.join(
                    centroids_lazy.select(["sequence_id"]), on="sequence_id", how="semi"
                ).select(new_df.columns)

                # Add the 'count' column to the filtered DataFrame (only for matching rows)
                if keep_cluster_sizes:
                    filtered_df = (
                        filtered_df.join(
                            centroids_lazy.select(["sequence_id", "count"]),
                            on="sequence_id",
                            how="left",
                        )
                        .with_columns(
                            pl.col("count").fill_null(
                                0
                            )  # Replace null values with 0 if necessary
                        )
                        .rename(
                            {
                                "count": "umi_count" if umi else "cluster_size"
                            }  # Rename the column
                        )
                    )

                # Sink the resulting DataFrame to a TSV file
                filtered_df.collect().write_csv(tsv_file, separator="\t")
                print(f"Output written to TSV file: {tsv_file}\n")

    if pool:
        heavies = pl.concat(pooled_heavies)
        lights = pl.concat(pooled_lights)

        all_consentroids = []

        # Processing pooled heavy chains
        heavy_consentroids = _process_chain_group(
            heavies,
            "Heavy",
            umi,
            min_cluster_size,
            clustering_threshold,
            consentroid,
            cluster_sizes_separator,
            keep_cluster_sizes,
            output_format,
            debug,
        )
        all_consentroids.extend(heavy_consentroids)

        # Processing pooled light chains
        light_consentroids = _process_chain_group(
            lights,
            "Lights",
            umi,
            min_cluster_size,
            clustering_threshold,
            consentroid,
            cluster_sizes_separator,
            keep_cluster_sizes,
            output_format,
            debug,
        )
        all_consentroids.extend(light_consentroids)

        if output_format == "fasta":
            # Saving pooled consentroids to FASTA file
            fasta_file = os.path.join(project_folder, "reduced_pool.fasta")
            to_fasta(all_consentroids, fasta_file)

        elif output_format == "airr":
            # Saving consentroids to TSV file
            tsv_file = os.path.join(project_folder, "reduced_pool.tsv")

            # Convert centroids to LazyFrame for more efficient processing
            centroids_lazy = pl.LazyFrame(all_consentroids)

            # Scan the complete original AIRR table
            new_df = pl.concat(
                [
                    pl.scan_csv(
                        f, separator="\t", low_memory=True if large_files else False
                    )
                    for f in files
                ]
            )

            # Filter new_df to only include rows with matching sequence_ids
            filtered_df = new_df.join(
                centroids_lazy.select(["sequence_id"]), on="sequence_id", how="semi"
            ).select(new_df.columns)

            # Add the 'count' column to the filtered DataFrame (only for matching rows)
            if keep_cluster_sizes:
                filtered_df = (
                    filtered_df.join(
                        centroids_lazy.select(["sequence_id", "count"]),
                        on="sequence_id",
                        how="left",
                    )
                    .with_columns(
                        pl.col("count").fill_null(
                            0
                        )  # Replace null values with 0 if necessary
                    )
                    .rename(
                        {
                            "count": "umi_count" if umi else "cluster_size"
                        }  # Rename the column
                    )
                )

            # Sink the resulting DataFrame to a TSV file
            filtered_df.collect().write_csv(tsv_file, separator="\t")
            print(f"Output written to TSV file: {tsv_file}\n")

    end = time.time()
    duration = end - start
    print(
        f"Data reduction complete. Total elapsed time: {datetime.timedelta(seconds=int(duration))}."
    )

    return
