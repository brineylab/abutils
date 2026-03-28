#!/usr/bin/env python
# filename: cluster.py


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
from typing import Iterable, Optional, Union

from ..bin import get_path as get_binary_path
from ..core.sequence import Sequence
from ..io import make_dir, to_fasta

__all__ = ["mmseqs_search"]


def mmseqs_search(
    query: Union[str, Sequence, Iterable[Sequence]],
    target: Union[str, Sequence, Iterable[Sequence]],
    output_path: str,
    temp_directory: Optional[str] = None,
    search_type: Optional[int] = None,
    max_seqs: Optional[int] = None,
    max_evalue: Optional[float] = None,
    sensitivity: Optional[float] = None,
    format_mode: Optional[int] = None,
    format_output: Optional[str] = None,
    verbosity: Optional[int] = None,
    threads: Optional[int] = None,
    split_memory_limit: Optional[str] = None,
    additional_cli_args: Optional[str] = None,
    mmseqs_bin: Optional[str] = None,
    log_to: Optional[str] = None,
    debug: bool = False,
) -> str:
    """
    Run MMseqs2_ easy-search.

    Parameters
    ----------
    query : Union[str, Sequence, Iterable[Sequence]]
        Query sequences. Can be either:
            - A path to a FASTA file
            - A path to an MMseqs2 database
            - A Sequence object
            - An iterable of Sequence objects

    target : Union[str, Sequence, Iterable[Sequence]]
        Target sequences. Can be either:
            - A path to a FASTA file
            - A path to an MMseqs2 database
            - A Sequence object
            - An iterable of Sequence objects

    output_path : str
        Path to save the output to. If the parent directory does not exist, it will be created.

    temp_directory : Optional[str]
        Path to save temporary files to. If not provided, the system temp directory will be used.

    search_type : int
        Search type. Options are:
          - 0: auto
          - 1: amino acid
          - 2: translated
          - 3: nucleotide
          - 4: translated nucleotide alignment
        If not provided, the default search type is 0 (auto).

    max_seqs : int
        Maximum number of matches to return for each query sequence. If not provided, the default is 300.

    max_evalue : float
        Only return matches below this E-value. If not provided, the default is 1.0e-3.

    sensitivity : float
        Sensitivity. 1.0 is fastest but least sensitive, 7.5 is slowest but most sensitive.
        If not provided, the default sensitivity is 5.7.

    format_mode : int
        Format mode. Options are:
          - 0: BLAST-TAB
          - 1: SAM
          - 2: BLAST-TAB + query/db length
          - 3: pretty HTML
          - 4: BLAST-TAB + column headers
        BLAST-TAB (0) and BLAST-TAB + column headers (4) support custom output formats using ``format_output``.
        If not provided, the default format mode is 0 (BLAST-TAB).

    format_output : str
        Comma separated list of output columns. Options are:
          - query
          - target
          - evalue
          - gapopen
          - pident
          - fident
          - nident
          - qstart
          - qend
          - qlen
          - tstart
          - tend
          - tlen
          - alnlen
          - raw
          - bits
          - cigar
          - qseq
          - tseq
          - qheader
          - theader
          - qaln
          - taln
          - qframe
          - mismatch
          - qcov
          - tcov
        Default is ``"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"``

    verbosity : int
        Verbosity level. 0 is quiet, 1 is +errors, 2 is +warnings, 3 is +info.
        If not provided, the default verbosity level is 3 (+info).

    threads : int
        Number of threads to use. If not provided, the default is to use all available threads.

    split_memory_limit : str
        Memory limit for split searches. Should be a string like ``"10G"`` or ``"100M"``.
        If not provided, the limit is the default set in ``mmseqs2`` (all system memory).

    additional_cli_args : str
        Additional command line arguments (as a string) which will be passed directly to ``mmseqs easy-search``.

        .. warning::
            Providing command-line options as both arguments to this function call and as `additional_cli_args`
            (for example, setting `max_seqs` to 100 and also providing ``"--max-seqs 100"`` to `additional_cli_args`)
            will cause an error. Don't do that.

    mmseqs_bin : str
        Path to the MMseqs binary. If not provided, the MMseqs binary packaged with ``abutils`` will be used.

    debug : bool
        If True, print the command and stdout/stderr to the console.

    Returns
    -------
    output_path : str
        Path to the output file.


    .. _MMseqs2:
        https://github.com/soedinglab/MMseqs2
    """
    # process input data
    if not isinstance(query, str):
        if isinstance(query, Sequence):
            query = [query]
        query = to_fasta(query)
    if not isinstance(target, str):
        if isinstance(target, Sequence):
            target = [target]
        target = to_fasta(target)

    # set up output and temp directories
    make_dir(os.path.dirname(output_path))
    if temp_directory is None:
        temp_directory = "/tmp"

    # build mmseqs command
    if mmseqs_bin is None:
        mmseqs_bin = get_binary_path("mmseqs")
    mmseqs_cmd = f"{mmseqs_bin} easy-search {query} {target} {output_path} {temp_directory} --search-type {search_type}"
    if sensitivity is not None:
        mmseqs_cmd += f" -s {sensitivity}"
    if max_seqs is not None:
        mmseqs_cmd += f" --max-seqs {max_seqs}"
    if max_evalue is not None:
        mmseqs_cmd += f" -e {max_evalue}"
    if format_mode is not None:
        mmseqs_cmd += f" --format-mode {format_mode}"
    if format_output is not None:
        mmseqs_cmd += f" --format-output {format_output}"
    if verbosity is not None:
        mmseqs_cmd += f" -v {verbosity}"
    if threads is not None:
        mmseqs_cmd += f" --threads {threads}"
    if split_memory_limit is not None:
        mmseqs_cmd += f" --split-memory-limit {split_memory_limit}"
    if additional_cli_args is not None:
        mmseqs_cmd += f" {additional_cli_args}"

    # run MMseqs
    p = sp.Popen(mmseqs_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode != 0:
        raise RuntimeError(f"MMseqs command failed: {stderr.decode('utf-8')}")
    if debug:
        print("STDOUT:\n", stdout.decode("utf-8"))
        print("\nSTDERR:\n", stderr.decode("utf-8"))
    if log_to is not None:
        log_to = os.path.abspath(log_to)
        make_dir(os.path.dirname(log_to))
        with open(log_to, "w") as f:
            f.write("STDOUT:\n-----\n")
            f.write(stdout.decode("utf-8"))
            f.write("\n\nSTDERR:\n-----\n")
            f.write(stderr.decode("utf-8"))
    return output_path
