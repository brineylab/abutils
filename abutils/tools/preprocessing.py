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


import os
import subprocess as sp
import tempfile
from typing import Iterable, Optional, Union

from natsort import natsorted
from tqdm.auto import tqdm

from ..bin import get_path as get_binary_path
from ..io import concatenate_files, delete_files, list_files, make_dir, rename_file


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
        return "_".join(self.basename.split("_")[:-4])

    @property
    def number(self):
        return self.basename.split(".")[0].split("_")[-1]

    @property
    def read(self):
        return self.basename.split("_")[-2]

    @property
    def lane(self):
        return self.basename.split("_")[-3]

    @property
    def sample(self):
        return self.basename.split("_")[-4]


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
        algo: str = "vsearch",
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
        verbose: bool = False,
        debug: bool = False,
    ) -> str:
        groups = self._group_by_lane()
        n_groups = len(groups)
        if verbose:
            print(f"{self.name}")
            print("-" * len(self.name))
            if n_groups > 1:
                groups = tqdm(groups, desc="  - merging lanes")
        # merge function
        if algo.lower() == "fastp":
            merge_func = merge_fastqs_fastp
        else:
            raise ValueError(f"Invalid merge algorithm: {algo}. Must be 'fastp'.")
        # merge files
        merged_files = []
        for group in groups:
            r1, r2 = natsorted(group, key=lambda x: x.read)
            merged_file = merge_func(
                r1.path,
                r2.path,
                os.path.join(merged_directory, f"{self.name}.{format.lower}"),
                name=self.name,
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
                output_format=format,
                debug=debug,
            )
            merged_files.append(merged_file)
        self.merged_file = os.path.join(merged_directory, f"{self.name}.{format.lower}")
        if len(merged_files) > 1:
            if verbose:
                print("  - concatenating merged files")
            concatenate_files(merged_files, self.merged_file)
            delete_files(merged_files)
        else:
            rename_file(merged_files, self.merged_file)
        return self.merged_file

    def _group_by_lane(self):
        lane_dict = {}
        for f in self.files:
            if f.lane not in lane_dict:
                lane_dict[f.lane] = []
            lane_dict[f.lane].append(f)
        return [l[0] for l in natsorted(lane_dict.items(), key=lambda x: x[0])]


def merge_fastqs(
    files: Union[str, Iterable],
    output_directory: str,
    log_directory: Optional[str] = None,
    output_format: str = "fastq",
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

    schema : str, optional
        Schema of the file names. Must be 'illumina'. Default is 'illumina'.

    algo : str, optional
        Algorithm to use for merging. Must be 'fastp'. Default is 'fastp'.

    binary_path : str, optional
        Path to the merge algorithm binary. If not provided, the binary packaged with abutils will be used.

    merge_args : str, optional
        Additional arguments (as a string) to pass to the merge function.

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
    else:
        raise ValueError(f"Invalid schema: {schema}. Must be 'illumina'.")
    # group files by sample
    file_pairs = group_fastq_pairs(files, verbose=verbose)
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
    groups = [
        MergeGroup(name, group_files, verbose=verbose)
        for name, group_files in group_dict.items()
    ]
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
        cmd += " --fastaout {merged_file}"
    elif output_format.lower() == "fastq":
        cmd += " --fastqout {merged_file}"
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
    merged_file: str,
    binary_path: Optional[str] = None,
    minimum_overlap: int = 30,
    allowed_mismatches: int = 5,
    allowed_mismatch_percent: float = 20.0,
    trim_adapters: bool = True,
    adapter_file: Optional[str] = None,
    quality_trim: bool = True,
    window_size: int = 4,
    quality_cutoff: int = 20,
    name: Optional[str] = None,
    log_directory: Optional[str] = None,
    additional_args: Optional[str] = None,
    debug: bool = False,
    **kwargs,
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

    binary_path : str, optional
        Path to a fastp binary. If not provided, the fastp binary packaged with abutils will be used.

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
    out_dir = os.path.dirname(merged_file)
    make_dir(out_dir)

    # get the vsearch binary
    if binary_path is None:
        binary_path = get_binary_path("fastp")

    # compile the vsearch command
    cmd = f"{binary_path} -i '{forward}' -I '{reverse}' --merge --merged_out '{merged_file}'"
    cmd += f" --overlap_length_require {minimum_overlap}"
    cmd += f" --overlap_diff_limit {allowed_mismatches}"
    cmd += f" --overlap_diff_percent_limit {allowed_mismatch_percent}"

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
        cmd += f" --cut_tail_window_size {window_size}"
        cmd += f" --cut_tail_mean_quality {quality_cutoff}"
    else:
        cmd += " --disable_quality_filtering"

    # log
    if log_directory is None:
        log_directory = tempfile.mkdtemp()
    if name is None:
        name = ".".join(os.path.basename(merged_file).split(".")[:-1])
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
    return merged_file
