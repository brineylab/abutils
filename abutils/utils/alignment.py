#!/usr/bin/env python
# filename: alignment.py
#
#
# Copyright (c) 2023 Bryan Briney
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


from abc import ABC, abstractmethod
from copy import copy, deepcopy
from io import StringIO
from itertools import groupby
import os
import platform
import subprocess as sp
import sys
import tempfile
from typing import Union, Iterable, Optional, Callable
import uuid

import parasail

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.SeqRecord import SeqRecord

from .decorators import lazy_property

from ..core.sequence import Sequence
from ..io import to_fasta


__all__ = [
    "mafft",
    "muscle",
    "muscle_v3",
    "local_alignment",
    "global_alignment",
    "semiglobal_alignment",
    "dot_alignment",
    "LocalAlignment",
    "GlobalAlignment",
    "SemiGlobalAlignment",
]


# -------------------------------------
#
#     MULTIPLE SEQUENCE ALIGNMENT
#
# -------------------------------------


def mafft(
    sequences: Union[str, Iterable],
    alignment_file: Optional[str] = None,
    fmt: str = "fasta",
    threads: int = -1,
    as_file: bool = False,
    as_string: bool = False,
    reorder: bool = True,
    mafft_bin: Optional[str] = None,
    id_key: Optional[str] = None,
    seq_key: Optional[str] = None,
    debug: bool = False,
    fasta: Optional[str] = None,
) -> Union[MultipleSeqAlignment, str]:
    """
    Performs multiple sequence alignment with `MAFFT`_.

    Parameters
    ----------
    sequences : (str, iterable)
        Can be one of several things:
            1. path to a FASTA-formatted file
            2. a FASTA-formatted string
            3. a list of BioPython ``SeqRecord`` objects
            4. a list of abutils ``Sequence`` objects
            5. a list of lists/tuples, of the format ``[sequence_id, sequence]``

    alignment_file : str, optional
        Path for the output alignment file. Required if ``as_file`` is ``True``.

    fmt : str, default='fasta'
        Format of the output alignment. Choices are 'fasta', 'phylip', and 'clustal'.
        Default is 'fasta'.

    threads : int, default=-1
        Number of threads for MAFFT to use. Default is ``-1``, which uses all
        available CPUs.

    as_file: bool, default=False
        If ``True``, returns the path to the alignment file. If ``False``,
        returns either a BioPython ``MultipleSeqAlignment`` object or the alignment
        output as a ``str``, depending on `as_string`. If `alignment_file` is not
        provided, a temporary file will be created with ``tempfile.NamedTemporaryFile``.
        Note that this temporary file is created in ``"/tmp"`` and may be removed
        by the operating system without notice.

    as_string: bool, default=False
        If ``True``, returns a the alignment output as a string. If ``False``,
        returns a BioPython ``MultipleSeqAlignment`` object (obtained by calling
        ``Bio.AlignIO.read()`` on the alignment file).

    mafft_bin : str, optional
        Path to a MAFFT executable. Default is ``None``, which results in
        using the default system MAFFT binary (just calling ``'mafft'``).

    id_key : str, default=None
        Key to retrieve the sequence ID. If not provided or missing, ``Sequence.id`` is used.

    sequence_key : str, default=None
        Key to retrieve the sequence. If not provided or missing, ``Sequence.sequence`` is used.

    debug : bool, default=False
        If ``True``, prints MAFFT's standard output and standard error.
        Default is ``False``.

    fasta : str, optional
        Path to a FASTA-formatted input file. Depricated (use `sequences`
        for all types if input), but retained for backwards compatibility.


    Returns
    -------
    alignment : str or ``MultipleSeqAlignment``
        If ``as_file`` is ``True``, returns a path to the output alignment file,
        Otherwise, returns a BioPython ``MultipleSeqAlignment`` object.


    .. _MAFFT:
        https://mafft.cbrc.jp/alignment/software/
    """
    # process input
    if fasta is not None:
        sequences = fasta
    ffile = to_fasta(sequences, id_key=id_key, sequence_key=seq_key)
    # configure output path
    if alignment_file is None:
        # as_file = False
        alignment_file = tempfile.NamedTemporaryFile(delete=False).name
    else:
        alignment_file = os.path.abspath(alignment_file)
    # do the alignment
    aln_format = ""
    if fmt.lower() == "clustal":
        aln_format = "--clustalout "
    if fmt.lower() == "phylip":
        aln_format = "--phylipout "
    if reorder:
        aln_format += "--reorder "
    if mafft_bin is None:
        mafft_bin = "mafft"
    mafft_cline = "{} --thread {} {}{} > {}".format(
        mafft_bin, threads, aln_format, ffile, alignment_file
    )
    mafft = sp.Popen(
        str(mafft_cline),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True,
        shell=True,
    )
    stdout, stderr = mafft.communicate()
    # return output
    if debug:
        print(mafft_cline)
        print(stdout)
        print(stderr)
    if os.stat(alignment_file).st_size == 0:
        err = "WARNING: output alignment file is empty. "
        err += "Verify that MAFFT is installed and the input data is valid,"
        print(err)
        return None
    if as_file:
        return alignment_file
    with open(alignment_file, "r") as f:
        if as_string:
            aln = f.read()
        else:
            aln = AlignIO.read(f, fmt)
    os.unlink(alignment_file)
    return aln


def muscle(
    sequences: Union[str, Iterable],
    alignment_file: Optional[str] = None,
    as_file: bool = False,
    as_string: bool = False,
    muscle_bin: Optional[str] = None,
    threads: Optional[int] = None,
    id_key: Optional[str] = None,
    seq_key: Optional[str] = None,
    debug: bool = False,
    fasta: Optional[str] = None,
) -> Union[MultipleSeqAlignment, str]:
    """
    Performs multiple sequence alignment with `MUSCLE`_ (version 5).

    Parameters
    ----------
    sequences : (str, iterable)
        Can be one of several things:
            1. path to a FASTA-formatted file
            2. a FASTA-formatted string
            3. a list of BioPython ``SeqRecord`` objects
            4. a list of abutils ``Sequence`` objects
            5. a list of lists/tuples, of the format ``[sequence_id, sequence]``

    alignment_file : str, optional
        Path for the output alignment file. Required if ``as_file`` is ``True``.

    as_file: bool, default=False
        If ``True``, returns the path to the alignment file. If ``False``,
        returns either a BioPython ``MultipleSeqAlignment`` object or the alignment
        output as a ``str``, depending on `as_string`. If `alignment_file` is not
        provided, a temporary file will be created with ``tempfile.NamedTemporaryFile``.
        Note that this temporary file is created in ``"/tmp"`` and may be removed
        by the operating system without notice.

    as_string: bool, default=False
        If ``True``, returns a the alignment output as a string. If ``False``,
        returns a BioPython ``MultipleSeqAlignment`` object (obtained by calling
        ``Bio.AlignIO.read()`` on the alignment file).

    muscle_bin : str, optional
        Path to a MUSCLE executable. If not provided, the MUSCLE binary bundled
        with ``abutils`` will be used.

    threads : int, default=None
        Number of threads for MUSCLE to use. If not provided, MUSCLE uses all
        available CPUs.

    id_key : str, default=None
        Key to retrieve the sequence ID. If not provided or missing, ``Sequence.id`` is used.

    sequence_key : str, default=None
        Key to retrieve the sequence. If not provided or missing, ``Sequence.sequence`` is used.

    debug : bool, default=False
        If ``True``, prints MAFFT's standard output and standard error.
        Default is ``False``.

    fasta : str, optional
        Path to a FASTA-formatted input file. Depricated (use `sequences`
        for all types if input), but retained for backwards compatibility.


    Returns
    -------
    alignment : str or ``MultipleSeqAlignment``
        If ``as_file`` is ``True``, returns a path to the output alignment file,
        Otherwise, returns a BioPython ``MultipleSeqAlignment`` object.


    .. _MUSCLE:
        https://www.drive5.com/muscle/
    """
    # process input
    if fasta is not None:
        sequences = fasta
    ffile = to_fasta(sequences, id_key=id_key, sequence_key=seq_key)
    # configure output path
    if alignment_file is None:
        # as_file = False
        alignment_file = tempfile.NamedTemporaryFile(delete=False).name
    else:
        alignment_file = os.path.abspath(alignment_file)
    # muscle binary
    if muscle_bin is None:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        muscle_bin = os.path.join(
            mod_dir, "bin/muscle_{}".format(platform.system().lower())
        )
    # do the alignment
    muscle_cline = f"{muscle_bin} -align {ffile} -output {alignment_file}"
    if threads is not None:
        muscle_cline += f" -threads {threads}"
    muscle = sp.Popen(str(muscle_cline), stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = muscle.communicate()
    if debug:
        print("muscle binary path:", muscle_bin)
        print("muscle command:", muscle_cline)
        print(stdout)
        print(stderr)
    # output
    if os.stat(alignment_file).st_size == 0:
        err = "WARNING: output alignment file is empty. "
        err += "Please verify that the input data is valid"
        print(err)
        return None
    if as_file:
        return alignment_file
    with open(alignment_file, "r") as f:
        if as_string:
            aln = f.read()
        else:
            aln = AlignIO.read(f, "fasta")
    os.unlink(alignment_file)
    return aln


def muscle_v3(
    sequences: Union[str, Iterable, None] = None,
    alignment_file: Optional[str] = None,
    fmt: str = "fasta",
    as_file: bool = False,
    as_string: bool = False,
    muscle_bin: Optional[str] = None,
    maxiters: Optional[int] = None,
    diags: Optional[int] = None,
    gap_open: Optional[int] = None,
    gap_extend: Optional[int] = None,
    id_key: Optional[str] = None,
    seq_key: Optional[str] = None,
    debug: bool = False,
    fasta: Optional[str] = None,
) -> Union[MultipleSeqAlignment, str]:
    """
    Performs multiple sequence alignment with `MUSCLE`_ (version 3).

    Parameters
    ----------
    sequences : (str, iterable)
        Can be one of several things:
            1. path to a FASTA-formatted file
            2. a FASTA-formatted string
            3. a list of BioPython ``SeqRecord`` objects
            4. a list of abutils ``Sequence`` objects
            5. a list of lists/tuples, of the format ``[sequence_id, sequence]``

    alignment_file : str, optional
        Path for the output alignment file. Required if ``as_file`` is ``True``.

    fmt : str, default='fasta'
        Format of the output alignment. Choices are 'fasta', 'phylip', and 'clustal'.
        Default is 'fasta'.

    threads : int, default=-1
        Number of threads for MUSCLE to use. Default is ``-1``, which uses all
        available CPUs.

    as_file: bool, default=False
        If ``True``, returns the path to the alignment file. If ``False``,
        returns either a BioPython ``MultipleSeqAlignment`` object or the alignment
        output as a ``str``, depending on `as_string`. Requires that
        `alignment_file` is also provided.

    as_string: bool, default=False
        If ``True``, returns a the alignment output as a string. If ``False``,
        returns a BioPython ``MultipleSeqAlignment`` object (obtained by calling
        ``Bio.AlignIO.read()`` on the alignment file).

    maxiters : int, default=-1
        Passed directly to MUSCLE using the ``-maxiters`` flag.

    diags : int, default=None
        Passed directly to MUSCLE using the ``-diags`` flag.

    gap_open : float, default=None
        Passed directly to MUSCLE using the ``-gapopen`` flag. Ignored
        if ``gap_extend`` is not also provided.

    gap_extend : float, default=None
        Passed directly to MUSCLE using the ``-gapextend`` flag. Ignored
        if ``gap_open`` is not also provided.

    muscle_bin : str, optional
        Path to a MUSCLE executable. If not provided, the MUSCLE binary bundled
        with ``abutils`` will be used.

    debug : bool, default=False
        If ``True``, prints MAFFT's standard output and standard error.
        Default is ``False``.

    fasta : str, optional
        Path to a FASTA-formatted input file. Depricated (use `sequences`
        for all types if input), but retained for backwards compatibility.


    Returns
    -------
    alignment : str or ``MultipleSeqAlignment``
        If ``as_file`` is ``True``, returns a path to the output alignment file,
        Otherwise, returns a BioPython ``MultipleSeqAlignment`` object.


    .. _MUSCLE:
        https://drive5.com/muscle/downloads_v3.htm

    """
    # process input
    if fasta is not None:
        sequences = fasta
    fasta_string = to_fasta(
        sequences,
        as_string=True,
        id_key=id_key,
        sequence_key=seq_key,
    )
    # configure output path
    if alignment_file is None:
        as_file = False
        alignment_file = tempfile.NamedTemporaryFile(delete=False).name
    else:
        alignment_file = os.path.abspath(alignment_file)
    # muscle binary
    if muscle_bin is None:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        muscle_bin = os.path.join(
            mod_dir, "bin/muscle3_{}".format(platform.system().lower())
        )
    # do the alignment
    aln_format = ""
    if fmt == "clustal":
        aln_format = " -clwstrict"
    muscle_cline = "{}{} ".format(muscle_bin, aln_format)
    if maxiters is not None:
        muscle_cline += " -maxiters {}".format(maxiters)
    if diags:
        muscle_cline += " -diags"
    if all([gap_open is not None, gap_extend is not None]):
        muscle_cline += " -gapopen {} -gapextend {}".format(gap_open, gap_extend)
    if debug:
        print("muscle binary path:", muscle_bin)
        print("muscle command:", muscle_cline)
    muscle = sp.Popen(
        str(muscle_cline),
        stdin=sp.PIPE,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True,
        shell=True,
    )
    alignment, stderr = muscle.communicate(input=fasta_string)
    if debug:
        print(stderr)
    # output
    if as_string:
        return alignment
    elif as_file:
        with open(alignment_file, "w") as f:
            f.write(alignment)
    else:
        return AlignIO.read(StringIO(alignment), fmt)


def consensus(aln, name=None, threshold=0.51, ambiguous="N"):
    summary_align = AlignInfo.SummaryInfo(aln)
    consensus = summary_align.gap_consensus(threshold=threshold, ambiguous=ambiguous)
    if name is None:
        name = uuid.uuid4()
    consensus_string = str(consensus).replace("-", "")
    return (name, consensus_string.upper())


# ----------------------------
#
#     PAIRWISE ALIGNMENT
#
# ----------------------------


class PairwiseAlignment(ABC):
    """
    Base class for local, global, and semi-global pairwise alignments.

    .. note::

        All comparisons between ``PairwiseAlignment``s are done on the
        ``score`` attribute. This was done so that sorting alignments
        like so::

            list_of_alignments.sort(reverse=True)

        produces a sorted list of alignments from the highest alignment
        score to the lowest.

    Attributes
    ----------

        query (Sequence): The input query sequence, as an AbTools
            ``Sequence`` object.

        target (Sequence): The input target sequence, as an AbTools
            ``Sequence`` object.

        target_id (str): ID of the target sequence.

        raw_query: The raw query, before conversion to a ``Sequence``.

        raw_target: The raw target, before conversion to a ``Sequence``.

    """

    def __init__(
        self,
        query: Sequence,
        target: Sequence,
        match: int = 3,
        mismatch: int = -2,
        gap_open: int = 5,
        gap_extend: int = 2,
        matrix: Optional[parasail.Matrix] = None,
        gap_open_penalty: Optional[int] = None,
        gap_extend_penalty: Optional[int] = None,
    ):
        self.query = query
        self.target = target
        self.match = int(match)
        self.mismatch = int(mismatch)
        self.gap_open = self.get_gap_open(gap_open, gap_open_penalty)
        self.gap_extend = self.get_gap_extend(gap_extend, gap_extend_penalty)
        self.matrix = self.get_matrix(matrix)
        self._alignment = None
        self._traceback = None
        self._cigar = None

        # for backwards compatibility:
        self.raw_query = query
        self.raw_target = target

    def __repr__(self):
        if len(self.aligned_query) > 20:
            qstring = "{}...{}".format(
                self.aligned_query[:10], self.aligned_query[-10:]
            )
            mstring = "{}...{}".format(
                self.alignment_midline[:10], self.alignment_midline[-10:]
            )
            tstring = "{}...{}".format(
                self.aligned_target[:10], self.aligned_target[-10:]
            )
        else:
            qstring = self.aligned_query
            mstring = self.alignment_midline
            tstring = self.aligned_target
        lines = []
        lines.append("Pairwise Alignment")
        lines.append("------------------\n")
        lines.append(f"query:  {qstring}")
        lines.append(f"        {mstring}")
        lines.append(f"target: {tstring}\n")
        lines.append(f"score: {self.score}")
        lines.append(f"type: {self.alignment_type}")
        lines.append(f"length: {len(self.aligned_query)}")
        print("\n".join(lines))
        return ""

    def __str__(self):
        lines = []
        lines.append(self.aligned_query)
        lines.append(self.alignment_midline)
        lines.append(self.aligned_target)
        return "\n".join(lines)

    def __len__(self):
        try:
            return len(self.aligned_query)
        except AttributeError:
            return 0

    def __eq__(self, other):
        if not hasattr(other, "score"):
            if type(other) in [int, float]:
                return self.score == other
            return False
        return self.score == other.score

    def __lt__(self, other):
        if not hasattr(other, "score"):
            if type(other) in [int, float]:
                return self.score < other
            return False
        return self.score < other.score

    def __le__(self, other):
        if not hasattr(other, "score"):
            if type(other) in [int, float]:
                return self.score <= other
            return False
        return self.score <= other.score

    def __gt__(self, other):
        if not hasattr(other, "score"):
            if type(other) in [int, float]:
                return self.score > other
            return False
        return self.score > other.score

    def __ge__(self, other):
        if not hasattr(other, "score"):
            if type(other) in [int, float]:
                return self.score >= other
            return False
        return self.score >= other.score

    @lazy_property
    def profile(self):
        return parasail.profile_create_16(
            self.query.sequence, self.match, self.mismatch, self.matrix
        )

    @lazy_property
    def alignment(self):
        return self.align()

    @property
    def cigar(self):
        cigar_string = self.alignment.cigar.decode.decode("utf-8")
        return CIGAR(cigar_string=cigar_string)

    @property
    def traceback(self):
        return self.alignment.traceback

    @lazy_property
    def score(self):
        return self.alignment.score

    @lazy_property
    def aligned_query(self):
        return self.traceback.query

    @lazy_property
    def aligned_target(self):
        return self.traceback.ref

    @lazy_property
    def alignment_midline(self):
        return self.traceback.comp.replace(".", " ")

    @lazy_property
    def query_begin(self):
        # return self.alignment.cigar.beg_query
        query_begin = self.alignment.cigar.beg_query
        if self.cigar[0].element == "I":
            query_begin += self.cigar[0].length
        return query_begin

    @lazy_property
    def query_end(self):
        return self.alignment.end_query

    @lazy_property
    def target_begin(self):
        # return self.alignment.cigar.beg_ref
        target_begin = self.alignment.cigar.beg_ref
        if self.cigar[0].element == "D":
            target_begin += self.cigar[0].length
        return target_begin

    @lazy_property
    def target_end(self):
        return self.alignment.end_ref

    @property
    def target_id(self):
        return self.target.id

    def align(self):
        """
        docstring for align
        """
        aln = self.alignment_function(
            self.query.sequence,
            self.target.sequence,
            self.gap_open,
            self.gap_extend,
            self.matrix,
        )
        return aln

    def get_gap_open(self, gap_open, gap_open_penalty):
        if gap_open_penalty is not None:
            gap_open = gap_open_penalty
        return abs(int(gap_open))

    def get_gap_extend(self, gap_extend, gap_extend_penalty):
        if gap_extend_penalty is not None:
            gap_extend = gap_extend_penalty
        return abs(int(gap_extend))

    def get_matrix(self, matrix):
        if matrix is None:
            residues = "ACDEFGHIKLMNPQRSTVWY"
            return parasail.matrix_create(residues, self.match, self.mismatch)
        elif isinstance(matrix, str):
            if os.path.isfile(matrix):
                try:
                    return parasail.Matrix(matrix)
                except:
                    err = (
                        f"\nERROR: the provided matrix file:\n{matrix}\ndoes not exist "
                    )
                    err += "or is improperly formatted. Please provide a valid matrix file.\n"
                    print(err)
                    sys.exit(1)
            else:
                try:
                    return parasail.Matrix(matrix.lower())
                except ValueError:
                    err = f"\nERROR: the provided matrix name {matrix} was not found. "
                    err += "Please provide a valid parasail matrix name.\n"
                    print(err)
                    sys.exit(1)
        elif isinstance(matrix, parasail.Matrix):
            return matrix
        else:
            err = "\nERROR: If supplied, matrix must be either the name of a built-in "
            err += "parasail matrix, the path to a properly formatted matrix file, "
            err += f"or a parasail.Matrix object. You supplied:\n{matrix}\n"
            print(err)
            sys.exit(1)


class LocalAlignment(PairwiseAlignment):
    """
    docstring for LocalAlignment
    """

    def __init__(
        self,
        query: Sequence,
        target: Sequence,
        match: int = 3,
        mismatch: int = -2,
        gap_open: int = 5,
        gap_extend: int = 2,
        matrix: Union[str, parasail.Matrix, None] = None,
        gap_open_penalty: Optional[int] = None,
        gap_extend_penalty: Optional[int] = None,
        alignment_function: Callable = parasail.sw_trace_striped_16,
    ):
        super().__init__(
            query=query,
            target=target,
            matrix=matrix,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
        )
        self.alignment_function = alignment_function
        self.alignment_type = "local"


class GlobalAlignment(PairwiseAlignment):
    """
    docstring for LocalAlignment
    """

    def __init__(
        self,
        query: Sequence,
        target: Sequence,
        match: int = 3,
        mismatch: int = -2,
        gap_open: int = 5,
        gap_extend: int = 2,
        matrix: Union[str, parasail.Matrix, None] = None,
        gap_open_penalty: Optional[int] = None,
        gap_extend_penalty: Optional[int] = None,
        alignment_function: Callable = parasail.nw_trace_striped_16,
    ):
        super().__init__(
            query=query,
            target=target,
            matrix=matrix,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
        )
        self.alignment_function = alignment_function
        self.alignment_type = "global"


class SemiGlobalAlignment(PairwiseAlignment):
    """
    docstring for LocalAlignment
    """

    def __init__(
        self,
        query: Sequence,
        target: Sequence,
        match: int = 3,
        mismatch: int = -2,
        gap_open: int = 5,
        gap_extend: int = 2,
        matrix: Union[str, parasail.Matrix, None] = None,
        gap_open_penalty: Optional[int] = None,
        gap_extend_penalty: Optional[int] = None,
        alignment_function: Callable = parasail.sg_trace_striped_16,
    ):
        super().__init__(
            query=query,
            target=target,
            matrix=matrix,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
        )
        self.alignment_function = alignment_function
        self.alignment_type = "semi-global"


class CIGAR:
    """
    Parsable representation of a CIGAR string. The CIGAR parsing code was taken
    almost vertabim from the `Cigar`_ package. Init argument is a CIGAR string.
    Iterating over a ``CIGAR`` object yields ``CIGARElement``s. The length of
    a ``CIGAR`` object is the length of the alignment encoded by the CIGAR string,
    not the length of the CIGAR string itself::

        cig = CIGAR("99=2X99=")
        len(cig) # --> 200
        len(cig.cigar_string) # --> 8

    .. _Cigar
        https://github.com/brentp/cigar
    """

    def __init__(self, cigar_string):
        self.cigar_string = cigar_string
        self._cigar_list = None

    def __str__(self):
        return self.cigar_string

    def __repr__(self):
        print(self.cigar_string)
        return ""

    def __len__(self):
        return sum([c.length for c in self.cigar_list])

    def __getitem__(self, slice):
        return self.cigar_list[slice]

    def __iter__(self):
        for c in self.cigar_list:
            yield c

    @property
    def cigar_list(self):
        if self._cigar_list is None:
            cigar_list = []
            cig_iter = groupby(self.cigar_string, lambda c: c.isdigit())
            for _, n in cig_iter:
                ce = CIGARElement(
                    length=int("".join(n)), element="".join(next(cig_iter)[1])
                )
                cigar_list.append(ce)
            self._cigar_list = cigar_list
        return self._cigar_list


class CIGARElement:
    """
    Representation of a single CIGAR element.
    """

    def __init__(self, length, element):
        self.length: int = int(length)
        self.element: str = element


def process_targets(target, targets):
    if target is not None:
        return [Sequence(target)]
    return [Sequence(t) for t in targets]


def local_alignment(
    query: Union[str, SeqRecord, Sequence, Iterable],
    target: Union[str, SeqRecord, Sequence, Iterable, None] = None,
    targets: Optional[Iterable] = None,
    match: int = 3,
    mismatch: int = -2,
    gap_open: int = -5,
    gap_extend: int = -2,
    matrix: Union[str, parasail.Matrix, None] = None,
    alignment_function: Callable = parasail.sw_trace_striped_16,
    aa: bool = False,
    gap_open_penalty: Optional[int] = None,
    gap_extend_penalty: Optional[int] = None,
) -> Union[LocalAlignment, Iterable]:
    """
    Striped Smith-Waterman local pairwise sequence alignment.

    Parameters
    ----------
    query : str, SeqRecord, Sequence, or list
        Query sequence. Can be any of the following:
            1. a nucleotide or amino acid sequence, as a string
            2. a Biopython ``SeqRecord`` object
            3. an abutils ``Sequence`` object
            4. a ``list`` or ``tuple`` of the format ``[seq_id, sequence]``

    target : str, SeqRecord, Sequence, or list, default=None
        Target sequence. Can be anything accepted by `query`. One of
        `target` or `targets` must be provided.

    target : str, SeqRecord, Sequence, or list, default=None
        Iterable of multiple target sequences. A ``list`` or ``tuple`` of
        anything accepted by `query`. One of `target` or `targets` must be provided.

    match : int, default=3
        Match score. Should be a positive integer.

    mismatch : int, default=-2
        Mismatch score. Typically should be less than or equal to ``0``.

    gap_open : int, default=-5
        Gap open score. Must be less than or equal to ``0``.

    gap_extend : int, default=-2
        Gap extension score. Must be less than or equal to ``0``.

    matrix : str, dict or parasail.Matrix, default=None
        Scoring matrix. Can be either:
            - the name of a ``parasail`` built-in matrix.
            - path to a matrix file (in a format accepted by ``parasail``).
            - a ``parasail.Matrix`` object
        If not provided, a matrix will be created using `match` and `mismatch`.

    alignment_function : Callable, default=parasail.sw_trace_striped_16
        ``parasail`` `local alignment function`_ to be used for alignment.

    aa : bool, default=False
        Not used. Retained for backwards compatibility.

    gap_open_penalty : int, default=None
        Depricated, but retained for backwards compatibility. Use `gap_open`
        insteaed.

    gap_extend_penalty : int, default=None
        Depricated, but retained for backwards compatibility. Use `gap_exted`
        insteaed.


    Returns
    -------
    alignment : LocalAlignment or iterable
        If a single target sequence is provided (via ``target``), a single ``LocalAlignment``
        object will be returned. If multiple target sequences are supplied (via ``targets``),
        a list of ``LocalAlignment`` objects will be returned.


    .. _local alignment function
        https://github.com/jeffdaily/parasail-python#standard-function-naming-convention
    """
    # input processing
    query = Sequence(query)
    targets = process_targets(target, targets)
    # alignment
    alignments = []
    for target in targets:
        aln = LocalAlignment(
            query=query,
            target=target,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix=matrix,
            alignment_function=alignment_function,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
        )
        alignments.append(aln)
    if len(alignments) == 1:
        return alignments[0]
    return alignments


def global_alignment(
    query: Union[str, SeqRecord, Sequence, Iterable],
    target: Union[str, SeqRecord, Sequence, Iterable, None] = None,
    targets: Optional[Iterable] = None,
    match: int = 3,
    mismatch: int = -2,
    gap_open: int = -5,
    gap_extend: int = -2,
    matrix: Union[str, parasail.Matrix, None] = None,
    alignment_function: Callable = parasail.nw_trace_striped_16,
    aa: bool = False,
    gap_open_penalty: Optional[int] = None,
    gap_extend_penalty: Optional[int] = None,
) -> Union[GlobalAlignment, Iterable]:
    """
    Needleman-Wunch global pairwise sequence alignment.

    Parameters
    ----------
    query : str, SeqRecord, Sequence, or list
        Query sequence. Can be any of the following:
            1. a nucleotide or amino acid sequence, as a string
            2. a Biopython ``SeqRecord`` object
            3. an abutils ``Sequence`` object
            4. a ``list`` or ``tuple`` of the format ``[seq_id, sequence]``

    target : str, SeqRecord, Sequence, or list, default=None
        Target sequence. Can be anything accepted by `query`. One of
        `target` or `targets` must be provided.

    target : str, SeqRecord, Sequence, or list, default=None
        Iterable of multiple target sequences. A ``list`` or ``tuple`` of
        anything accepted by `query`. One of `target` or `targets` must be provided.

    match : int, default=3
        Match score. Should be a positive integer.

    mismatch : int, default=-2
        Mismatch score. Typically should be less than or equal to ``0``.

    gap_open : int, default=-5
        Gap open score. Must be  less than or equal to ``0``.

    gap_extend : int, default=-2
        Gap extension score. Must be less than or equal to ``0``.

    matrix : str, dict or parasail.Matrix, default=None
        Scoring matrix. Can be either:
            - the name of a ``parasail`` built-in matrix.
            - path to a matrix file (in a format accepted by ``parasail``).
            - a ``parasail.Matrix`` object
        If not provided, a matrix will be created using `match` and `mismatch`.

    alignment_function : Callable, default=parasail.sw_trace_striped_16
        ``parasail`` `global alignment function`_ to be used for alignment.

    aa : bool, default=False
        Not used. Retained for backwards compatibility.

    gap_open_penalty : int, default=None
        Depricated, but retained for backwards compatibility. Use `gap_open`
        insteaed.

    gap_extend_penalty : int, default=None
        Depricated, but retained for backwards compatibility. Use `gap_exted`
        insteaed.


    Returns
    -------
    alignment : GlobalAlignment or iterable
        If a single target sequence is provided (via ``target``), a single ``GlobalAlignment``
        object will be returned. If multiple target sequences are supplied (via ``targets``),
        a list of ``GlobalAlignment`` objects will be returned.


    .. _global alignment function
        https://github.com/jeffdaily/parasail-python#standard-function-naming-convention
    """
    # input processing
    query = Sequence(query)
    targets = process_targets(target, targets)
    # alignment
    alignments = []
    for target in targets:
        aln = GlobalAlignment(
            query=query,
            target=target,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix=matrix,
            alignment_function=alignment_function,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
        )
        alignments.append(aln)
    if len(alignments) == 1:
        return alignments[0]
    return alignments


def semiglobal_alignment(
    query: Union[str, SeqRecord, Sequence, Iterable],
    target: Union[str, SeqRecord, Sequence, Iterable, None] = None,
    targets: Optional[Iterable] = None,
    match: int = 3,
    mismatch: int = -2,
    gap_open: int = -5,
    gap_extend: int = -2,
    matrix: Union[str, parasail.Matrix, None] = None,
    alignment_function: Callable = parasail.sg_trace_striped_16,
    aa: bool = False,
    gap_open_penalty: Optional[int] = None,
    gap_extend_penalty: Optional[int] = None,
) -> Union[SemiGlobalAlignment, Iterable]:
    """
    Semi-global pairwise sequence alignment.

    Parameters
    ----------
    query : str, SeqRecord, Sequence, or list
        Query sequence. Can be any of the following:
            1. a nucleotide or amino acid sequence, as a string
            2. a Biopython ``SeqRecord`` object
            3. an abutils ``Sequence`` object
            4. a ``list`` or ``tuple`` of the format ``[seq_id, sequence]``

    target : str, SeqRecord, Sequence, or list, default=None
        Target sequence. Can be anything accepted by `query`. One of
        `target` or `targets` must be provided.

    target : str, SeqRecord, Sequence, or list, default=None
        Iterable of multiple target sequences. A ``list`` or ``tuple`` of
        anything accepted by `query`. One of `target` or `targets` must be provided.

    match : int, default=3
        Match score. Should be a positive integer.

    mismatch : int, default=-2
        Mismatch score. Typically should be less than or equal to ``0``.

    gap_open : int, default=-5
        Gap open score. Must be  less than or equal to ``0``.

    gap_extend : int, default=-2
        Gap extension score. Must be less than or equal to ``0``.

    matrix : str, dict or parasail.Matrix, default=None
        Scoring matrix. Can be either:
            - the name of a ``parasail`` built-in matrix.
            - path to a matrix file (in a format accepted by ``parasail``).
            - a ``parasail.Matrix`` object
        If not provided, a matrix will be created using `match` and `mismatch`.

    alignment_function : Callable, default=parasail.sw_trace_striped_16
        ``parasail`` `semi-global alignment function`_ to be used for alignment.

    aa : bool, default=False
        Not used. Retained for backwards compatibility.

    gap_open_penalty : int, default=None
        Depricated, but retained for backwards compatibility. Use `gap_open`
        insteaed.

    gap_extend_penalty : int, default=None
        Depricated, but retained for backwards compatibility. Use `gap_exted`
        insteaed.


    Returns
    -------
    alignment : SemiGlobalAlignment or iterable
        If a single target sequence is provided (via ``target``), a single ``SemiGlobalAlignment``
        object will be returned. If multiple target sequences are supplied (via ``targets``),
        a list of ``SemiGlobalAlignment`` objects will be returned.


    .. _semi-global alignment function
        https://github.com/jeffdaily/parasail-python#standard-function-naming-convention
    """
    # input processing
    query = Sequence(query)
    targets = process_targets(target, targets)
    # alignment
    alignments = []
    for target in targets:
        aln = SemiGlobalAlignment(
            query=query,
            target=target,
            match=match,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            matrix=matrix,
            alignment_function=alignment_function,
            gap_open_penalty=gap_open_penalty,
            gap_extend_penalty=gap_extend_penalty,
        )
        alignments.append(aln)
    if len(alignments) == 1:
        return alignments[0]
    return alignments


# def local_alignment(
#     query: Union[str, SeqRecord, Sequence, Iterable],
#     target: Union[str, SeqRecord, Sequence, Iterable, None] = None,
#     targets: Optional[Iterable] = None,
#     match: int = 3,
#     mismatch: int = -2,
#     gap_open: int = -5,
#     gap_extend: int = -2,
#     matrix: Union[str, dict, None] = None,
#     aa: bool = False,
#     gap_open_penalty: Optional[int] = None,
#     gap_extend_penalty: Optional[int] = None,
# ):
#     """
#     Striped Smith-Waterman local pairwise alignment.

#     Parameters
#     ----------
#     query : str, SeqRecord, Sequence, or list
#         Query sequence. Can be any of the following:
#             1. a nucleotide or amino acid sequence, as a string
#             2. a Biopython ``SeqRecord`` object
#             3. an abutils ``Sequence`` object
#             4. a ``list`` or ``tuple`` of the format ``[seq_id, sequence]``

#     target : str, SeqRecord, Sequence, or list, default=None
#         Target sequence. Can be anything accepted by `query`. One of
#         `target` or `targets` must be provided.

#     target : str, SeqRecord, Sequence, or list, default=None
#         Iterable of multiple target sequences. A ``list`` or ``tuple`` of
#         anything accepted by `query`. One of `target` or `targets` must be provided.

#     match : int, default=3
#         Match score. Should be a positive integer.

#     mismatch : int, default=-2
#         Mismatch score. Typically should be less than or equal to ``0``.

#     gap_open : int, default=-5
#         Gap open score. Typically should be less than or equal to ``0``.

#     gap_extend : int, default=-2
#         Gap extension score. Typically should be less than or equal to ``0``.

#     matrix : str, dict, default=None
#         Scoring matrix for amino acid alignments. Can be either:
#             - The name of a built-in matrix. Current options are
#             ``"blosum62"`` and ``"pam250"``.
#             - A nested dictionary, giving an alignment score for each residue
#             pair. Should be formatted such that retrieving the alignment score
#             for ``"A"`` and ``"G"`` is accomplished by::
#                 matrix['A']['G']
#         If not provided, ``"blosum62"`` is used.

#     aa : bool, default=False
#         Must be set to ``True`` if sequences are amino acids.

#     gap_open_penalty : int, default=None
#         Depricated, but retained for backwards compatibility. Use `gap_open`
#         insteaed. Should be a positive integer.

#     gap_extend_penalty : int, default=None
#         Depricated, but retained for backwards compatibility. Use `gap_exted`
#         insteaed. Should be a positive integer.


#     Returns
#     -------
#     alignment : SSWAlignment or iterable
#         If a single target sequence is provided (via ``target``), a single ``SSWAlignment``
#         object will be returned. If multiple target sequences are supplied (via ``targets``),
#         a list of ``SSWAlignment`` objects will be returned.
#     """
#     if aa and matrix is None:
#         matrix = "blosum62"
#     if not any([target, targets]):
#         err = "\nERROR: You must supply either ``target`` or ``targets``.\n"
#         raise RuntimeError(err)
#     if target:
#         targets = [
#             target,
#         ]
#     # to maintain backward compatibility with earlier abutils API
#     if gap_open_penalty is not None:
#         gap_open = -gap_open_penalty
#     if gap_extend_penalty is not None:
#         gap_extend = -gap_extend_penalty
#     # do the alignment
#     alignments = []
#     for t in targets:
#         try:
#             alignment = SSWAlignment(
#                 query=query,
#                 target=t,
#                 match=match,
#                 mismatch=mismatch,
#                 matrix=matrix,
#                 gap_open=-gap_open,
#                 gap_extend=-gap_extend,
#                 aa=aa,
#             )
#             alignments.append(alignment)
#         except IndexError:
#             continue
#     if len(alignments) == 1:
#         return alignments[0]
#     return alignments


# def local_alignment_biopython(
#     query,
#     target=None,
#     targets=None,
#     match=3,
#     mismatch=-2,
#     matrix=None,
#     gap_open=-5,
#     gap_extend=-2,
#     aa=False,
# ):
#     if not target and not targets:
#         err = "ERROR: You must supply a target sequence (or sequences)."
#         raise RuntimeError(err)
#     if target:
#         targets = [
#             target,
#         ]
#     alignments = []
#     for t in targets:
#         try:
#             alignment = alignment = BiopythonAlignment(
#                 query=query,
#                 target=t,
#                 match=match,
#                 mismatch=mismatch,
#                 matrix=matrix,
#                 gap_open=gap_open,
#                 gap_extend=gap_extend,
#                 aa=aa,
#             )
#             alignments.append(alignment)
#         except IndexError:
#             continue
#     if len(alignments) == 1:
#         return alignments[0]
#     return alignments


# def global_alignment(
#     query,
#     target=None,
#     targets=None,
#     match=3,
#     mismatch=-2,
#     gap_open=-5,
#     gap_extend=-2,
#     score_match=None,
#     score_mismatch=None,
#     score_gap_open=None,
#     score_gap_extend=None,
#     matrix=None,
#     aa=False,
# ):
#     """
#     Needleman-Wunch global pairwise alignment.

#     With ``global_alignment``, you can score an alignment using different
#     paramaters than were used to compute the alignment. This allows you to
#     compute pure identity scores (match=1, mismatch=0) on pairs of sequences
#     for which those alignment parameters would be unsuitable. For example::

#         seq1 = 'ATGCAGC'
#         seq2 = 'ATCAAGC'

#     using identity scoring params (match=1, all penalties are 0) for both alignment
#     and scoring produces the following alignment::

#         ATGCA-GC
#         || || ||
#         AT-CAAGC

#     with an alignment score of 6 and an alignment length of 8 (identity = 75%). But
#     what if we want to calculate the identity of a gapless alignment? Using::

#         global_alignment(seq1, seq2,
#                          gap_open=20,
#                          score_match=1,
#                          score_mismatch=0,
#                          score_gap_open=10,
#                          score_gap_extend=1)

#     we get the following alignment::

#         ATGCAGC
#         ||  |||
#         ATCAAGC

#     which has an score of 5 and an alignment length of 7 (identity = 71%). Obviously,
#     this is an overly simple example (it would be much easier to force gapless alignment
#     by just iterating over each sequence and counting the matches), but there are several
#     real-life cases in which different alignment and scoring paramaters are desirable.

#     Args:

#         query: Query sequence. ``query`` can be one of four things:

#             1. a nucleotide or amino acid sequence, as a string

#             2. a Biopython ``SeqRecord`` object

#             3. an AbTools ``Sequence`` object

#             4. a list/tuple of the format ``[seq_id, sequence]``

#         target: A single target sequence. ``target`` can be anything that
#             ``query`` accepts.

#         targets (list): A list of target sequences, to be proccssed iteratively.
#             Each element in the ``targets`` list can be anything accepted by
#             ``query``.

#         match (int): Match score for alignment. Should be a positive integer. Default is 3.

#         mismatch (int): Mismatch score for alignment. Should be a negative integer. Default is -2.

#         gap_open (int): Penalty for opening gaps in alignment. Should be a negative integer.
#             Default is -5.

#         gap_extend (int): Penalty for extending gaps in alignment. Should be a negative
#             integer. Default is -2.

#         score_match (int): Match score for scoring the alignment. Should be a positive integer.
#             Default is to use the score from ``match`` or ``matrix``, whichever is provided.

#         score_mismatch (int): Mismatch score for scoring the alignment. Should be a negative
#             integer. Default is to use the score from ``mismatch`` or ``matrix``, whichever
#             is provided.

#         score_gap_open (int): Gap open penalty for scoring the alignment. Should be a negative
#             integer. Default is to use ``gap_open``.

#         score_gap_extend (int): Gap extend penalty for scoring the alignment. Should be a negative
#             integer. Default is to use ``gap_extend``.

#         matrix (str, dict): Alignment scoring matrix. Two options for passing the alignment matrix:

#             - The name of a built-in matrix. Current options are ``blosum62`` and ``pam250``.

#             - A nested dictionary, giving an alignment score for each residue pair. Should be
#               formatted such that retrieving the alignment score for A and G is accomplished by::

#                 matrix['A']['G']

#         aa (bool): Must be set to ``True`` if aligning amino acid sequences. Default
#             is ``False``.

#     Returns:

#         If a single target sequence is provided (via ``target``), a single ``NWAlignment``
#         object will be returned. If multiple target sequences are supplied (via ``targets``),
#         a list of ``NWAlignment`` objects will be returned.
#     """
#     if not target and not targets:
#         err = "ERROR: You must supply a target sequence (or sequences)."
#         raise RuntimeError(err)
#     if target:
#         targets = [
#             target,
#         ]
#     if type(targets) not in (list, tuple):
#         err = "ERROR: ::targets:: requires an iterable (list or tuple)."
#         err += "For a single sequence, use ::target::"
#         raise RuntimeError(err)
#     alignments = []
#     for t in targets:
#         alignment = NWAlignment(
#             query=query,
#             target=t,
#             match=match,
#             mismatch=mismatch,
#             gap_open=gap_open,
#             gap_extend=gap_extend,
#             score_match=score_match,
#             score_mismatch=score_mismatch,
#             score_gap_open=score_gap_open,
#             score_gap_extend=score_gap_extend,
#             matrix=matrix,
#             aa=aa,
#         )
#         alignments.append(alignment)
#     if target is not None:
#         return alignments[0]
#     return alignments


# class BaseAlignment(object):
#     """
#     Base class for local and global pairwise alignments.

#     .. note::

#         All comparisons between ``BaseAlignments``
#         are done on the ``score`` attribute (which must be implemented
#         by any classes that subclass ``BaseAlignment``). This was done
#         so that sorting alignments like so::

#             alignments = [list of alignments]
#             alignments.sort(reverse=True)

#         results in a sorted list of alignments from the highest alignment
#         score to the lowest.

#     Attributes:

#         query (Sequence): The input query sequence, as an AbTools
#             ``Sequence`` object.

#         target (Sequence): The input target sequence, as an AbTools
#             ``Sequence`` object.

#         target_id (str): ID of the target sequence.

#         raw_query: The raw query, before conversion to a ``Sequence``.

#         raw_target: The raw target, before conversion to a ``Sequence``.

#     """

#     def __init__(
#         self, query, target, matrix, match, mismatch, gap_open, gap_extend, aa
#     ):
#         super(BaseAlignment, self).__init__()
#         self.raw_query = query
#         self.raw_target = target
#         self.query = self._process_sequence(query, aa=aa)
#         self.target = self._process_sequence(target, aa=aa)
#         self._matrix = matrix
#         self._match = int(match)
#         self._mismatch = int(mismatch)
#         self._gap_open = int(gap_open)
#         self._gap_extend = int(gap_extend)
#         self._aa = bool(aa)

#     def __repr__(self):
#         if len(self.aligned_query) > 20:
#             qstring = "{}...{}".format(
#                 self.aligned_query[:10], self.aligned_query[-10:]
#             )
#             mstring = "{}...{}".format(
#                 self.alignment_midline[:10], self.alignment_midline[-10:]
#             )
#             tstring = "{}...{}".format(
#                 self.aligned_target[:10], self.aligned_target[-10:]
#             )
#         else:
#             qstring = self.aligned_query
#             mstring = self.alignment_midline
#             tstring = self.aligned_target
#         return_string = "\n\n"
#         return_string += "Pairwise Alignment\n"
#         return_string += "------------------\n\n"
#         return_string += "query:  {}\n".format(qstring)
#         return_string += "        {}\n".format(mstring)
#         return_string += "target: {}\n\n".format(tstring)
#         return_string += "score: {}\n".format(str(self.score))
#         return_string += "type: {}\n".format(self.alignment_type)
#         return_string += "length: {}".format(str(len(self.aligned_query)))
#         print(return_string)
#         return ""

#     def __str__(self):
#         return_string = ""
#         return_string += "{}\n".format(self.aligned_query)
#         return_string += "{}\n".format(self.alignment_midline)
#         return_string += "{}\n".format(self.aligned_target)
#         return return_string

#     def __len__(self):
#         return len(self.aligned_query)

#     def __eq__(self, other):
#         if not hasattr(other, "score"):
#             if type(other) in [int, float]:
#                 return self.score == other
#             return False
#         return self.score == other.score

#     def __lt__(self, other):
#         if not hasattr(other, "score"):
#             if type(other) in [int, float]:
#                 return self.score == other
#             return False
#         return self.score < other.score

#     def __le__(self, other):
#         if not hasattr(other, "score"):
#             if type(other) in [int, float]:
#                 return self.score == other
#             return False
#         return self.score <= other.score

#     def __gt__(self, other):
#         if not hasattr(other, "score"):
#             if type(other) in [int, float]:
#                 return self.score == other
#             return False
#         return self.score > other.score

#     def __ge__(self, other):
#         if not hasattr(other, "score"):
#             if type(other) in [int, float]:
#                 return self.score == other
#             return False
#         return self.score >= other.score

#     @property
#     def target_id(self):
#         return self._target_id

#     @target_id.setter
#     def target_id(self, target_id):
#         self._target_id = target_id

#     @property
#     def alignment_midline(self):
#         if self._alignment_midline is None:
#             self._alignment_midline = self._get_alignment_midline()
#         return self._alignment_midline

#     @staticmethod
#     def _process_sequence(sequence, aa):
#         if type(sequence) == Sequence:
#             return sequence
#         return Sequence(sequence)

#     def _get_alignment_midline(self):
#         midline = ""
#         for q, t in zip(self.aligned_query, self.aligned_target):
#             if q == t:
#                 midline += "|"
#             else:
#                 midline += " "
#         return midline


# class SSWAlignment(BaseAlignment):
#     """
#     Structure for performing and analyzing a Smith-Waterman local alignment.

#     .. note:

#         Exposed attributes and methods are the same as ``NWAlignment``, so
#         local and global alignmnts can be handled in the same way. In fact,
#         since comparisons are made based on score, local and global alignments
#         can be directly compared with constructions like::

#             local_aln == global_aln
#             local_aln > global_aln
#             alignments = sorted([global_aln, local_aln])

#     Attributes:

#         alignment_type (str): Is 'local' for all ``SSWAlignment`` objects.

#         aligned_query (str): The aligned query sequence (including gaps).

#         aligned_target (str): The aligned target sequence (including gaps).

#         alignment_midline (str): Midline for the aligned sequences, with ``|`` indicating
#           matches and a gap indicating mismatches::

#             print(aln.aligned_query)
#             print(aln.alignment_midline)
#             print(aln.aligned_target)

#             # ATGC
#             # || |
#             # ATCC

#         score (int): Alignment score.

#         query_begin (int): Position in the raw query sequence at which
#             the optimal alignment begins.

#         query_end (int): Position in the raw query sequence at which the
#             optimal alignment ends.

#         target_begin (int): Position in the raw target sequence at which
#             the optimal alignment begins.

#         target_end (int): Position in the raw target sequence at which the
#             optimal alignment ends.
#     """

#     def __init__(
#         self,
#         query,
#         target,
#         match=3,
#         mismatch=-2,
#         matrix=None,
#         gap_open=5,
#         gap_extend=2,
#         aa=False,
#     ):
#         super(SSWAlignment, self).__init__(
#             query, target, matrix, match, mismatch, gap_open, gap_extend, aa
#         )

#         self.alignment_type = "local"
#         self._alignment = self._align()
#         self._aligned_query = self._alignment.aligned_query_sequence
#         self._aligned_target = self._alignment.aligned_target_sequence
#         self._alignment_midline = None
#         self.score = self._alignment.optimal_alignment_score
#         self.cigar = self._alignment.cigar
#         self.query_begin = self._alignment.query_begin
#         self.query_end = self._alignment.query_end
#         self.target_begin = self._alignment.target_begin
#         self.target_end = self._alignment.target_end_optimal
#         self._alignment = None

#     @property
#     def alignment(self):
#         return self._alignment

#     @property
#     def aligned_query(self):
#         return self._aligned_query

#     @aligned_query.setter
#     def aligned_query(self, aligned_query):
#         self._aligned_query = aligned_query
#         self._alignment_midline = None

#     @property
#     def aligned_target(self):
#         return self._aligned_target

#     @aligned_target.setter
#     def aligned_target(self, aligned_target):
#         self._aligned_target = aligned_target
#         self._alignment_midline = None

#     def _align(self):
#         aligner = StripedSmithWaterman(
#             self.query.sequence,
#             match_score=self._match,
#             mismatch_score=self._mismatch,
#             gap_open_penalty=self._gap_open,
#             gap_extend_penalty=self._gap_extend,
#             substitution_matrix=self._matrix,
#             protein=self._aa,
#         )
#         return aligner(self.target.sequence)


# class BiopythonAlignment(BaseAlignment):
#     def __init__(
#         self,
#         query,
#         target,
#         match=3,
#         mismatch=-2,
#         matrix=None,
#         gap_open=5,
#         gap_extend=2,
#         aa=False,
#     ):
#         super(BiopythonAlignment, self).__init__(
#             query, target, matrix, match, mismatch, gap_open, gap_extend, aa
#         )

#         self.alignment_type = "local"
#         self._aln = self._align()
#         aln_query, aln_target, score, begin, end = self._aln
#         self._aligned_query = aln_query[begin:end]
#         self._aligned_target = aln_target[begin:end]
#         self._alignment_midline = None
#         self.score = score
#         self.query_begin = self._get_begin_pos(aln_query, begin)
#         self.query_end = self._get_end_pos(aln_query, end)
#         self.target_begin = self._get_begin_pos(aln_target, begin)
#         self.target_end = self._get_end_pos(aln_target, end)

#     @property
#     def alignment(self):
#         return self._aln

#     @property
#     def aligned_query(self):
#         return self._aligned_query

#     @aligned_query.setter
#     def aligned_query(self, aligned_query):
#         self._aligned_query = aligned_query
#         self._alignment_midline = None

#     @property
#     def aligned_target(self):
#         return self._aligned_target

#     @aligned_target.setter
#     def aligned_target(self, aligned_target):
#         self._aligned_target = aligned_target
#         self._alignment_midline = None

#     def _align(self):
#         aln = pairwise2.align.localms(
#             self.query.sequence,
#             self.target.sequence,
#             self._match,
#             self._mismatch,
#             self._gap_open,
#             self._gap_extend,
#         )
#         return aln[0]

#     def _get_begin_pos(self, seq, begin):
#         dashes = seq.count("-", 0, begin)
#         return begin - dashes

#     def _get_end_pos(self, seq, end):
#         return len(seq[:end].replace("-", ""))


# class NWAlignment(BaseAlignment):
#     """
#     Structure for performing and analyzing a Needleman-Wunch global alignment.

#     .. note:

#         Exposed attributes and methods are the same as ``SSWAlignment``, so
#         local and global alignmnts can be handled in the same way. In fact,
#         since comparisons are made based on score, local and global alignments
#         can be directly compared with constructions like::

#             local_aln == global_aln
#             local_aln > global_aln
#             alignments = sorted([global_aln, local_aln])

#     Attributes:

#         alignment_type (str): Is 'global' for all ``NWAlignment`` objects.

#         aligned_query (str): The aligned query sequence (including gaps).

#         aligned_target (str): The aligned target sequence (including gaps).

#         alignment_midline (str): Midline for the aligned sequences, with
#           ``|`` indicating matches and a gap indicating mismatches::

#             print(aln.aligned_query)
#             print(aln.alignment_midline)
#             print(aln.aligned_target)

#             # ATGC
#             # || |
#             # ATCC

#         score (int): Alignment score.

#         query_begin (int): Position in the raw query sequence at which
#             the optimal alignment begins.

#         query_end (int): Position in the raw query sequence at which the
#             optimal alignment ends.

#         target_begin (int): Position in the raw target sequence at which
#             the optimal alignment begins.

#         target_end (int): Position in the raw target sequence at which the
#             optimal alignment ends.
#     """

#     def __init__(
#         self,
#         query,
#         target,
#         match=3,
#         mismatch=-2,
#         gap_open=-5,
#         gap_extend=-2,
#         score_match=None,
#         score_mismatch=None,
#         score_gap_open=None,
#         score_gap_extend=None,
#         matrix=None,
#         aa=False,
#     ):
#         super(NWAlignment, self).__init__(
#             query, target, matrix, match, mismatch, gap_open, gap_extend, aa
#         )
#         self.alignment_type = "global"
#         self._score_match = int(score_match) if score_match is not None else None
#         self._score_mismatch = (
#             int(score_mismatch) if score_mismatch is not None else None
#         )
#         self._score_gap_open = (
#             int(score_gap_open) if score_gap_open is not None else None
#         )
#         self._score_gap_extend = (
#             int(score_gap_extend) if score_gap_extend is not None else None
#         )
#         self._matrix = matrix
#         self._alignment = self._align()
#         self._aligned_query = self._alignment[0]
#         self._aligned_target = self._alignment[1]
#         self._alignment_midline = None
#         self.score = self._score_alignment()

#     @property
#     def alignment(self):
#         return self._aln

#     @property
#     def aligned_query(self):
#         return self._aligned_query

#     @aligned_query.setter
#     def aligned_query(self, aligned_query):
#         self._aligned_query = aligned_query
#         self._alignment_midline = None

#     @property
#     def aligned_target(self):
#         return self._aligned_target

#     @aligned_target.setter
#     def aligned_target(self, aligned_target):
#         self._aligned_target = aligned_target
#         self._alignment_midline = None

#     def _get_matrix_file(self, match=None, mismatch=None, matrix=None):
#         matrix_dir = os.path.join(
#             os.path.dirname(os.path.abspath(__file__)), "matrices"
#         )
#         builtins = [os.path.basename(f) for f in list_files(matrix_dir)]
#         if matrix is None:
#             matrix_name = "match{}mismatch{}".format(abs(match), abs(mismatch))
#             if matrix_name not in builtins:
#                 self._build_matrix_from_params(
#                     match, mismatch, os.path.join(matrix_dir, matrix_name)
#                 )
#             return os.path.join(matrix_dir, matrix_name)
#         elif isinstance(matrix, dict):
#             dhash = hashlib.md5()
#             encoded = json.dumps(matrix, sort_keys=True).encode()
#             dhash.update(encoded)
#             matrix_name = dhash.hexdigest()
#             if matrix_name not in builtins:
#                 self._build_matrix_from_dict(
#                     matrix, os.path.join(matrix_dir, matrix_name)
#                 )
#             return os.path.join(matrix_dir, matrix_name)
#         else:
#             if self._matrix.lower() in builtin_names:
#                 return os.path.join(matrix_dir, self._matrix.lower())
#             else:
#                 err = "The supplied matrix name ({}) does not exist. ".format(matrix)
#                 err += "Built-in matrices are: {}".format("\n".join(builtins))
#                 raise RuntimeError(err)

#         # if matrix is not None:
#         #     matrix_name = matrix
#         # else:
#         #     matrix_name = 'match{}mismatch{}'.format(abs(match), abs(mismatch))
#         # if matrix_name.lower() in builtins:
#         #     return os.path.join(matrix_dir, matrix_name)
#         # builtin_names = [os.path.basename(f) for f in list_files(matrix_dir)]
#         # if self._matrix is not None:
#         #     if self._matrix.lower() in builtin_names:
#         #         return os.path.join(matrix_dir, self._matrix.lower())
#         #     else:
#         #         err = 'The supplied matrix name ({}) does not exist. '.format(matrix)
#         #         err += 'Built-in matrices are: {}'.format(', '.join(builtins))
#         #         raise RuntimeError(err)
#         # else:
#         #     self._build_matrix_from_params(match, mismatch, os.path.join(matrix_dir, matrix_name))
#         #     return os.path.join(matrix_dir, matrix_name)

#     def _align(self):
#         matrix = self._get_matrix_file(
#             match=self._match, mismatch=self._mismatch, matrix=self._matrix
#         )
#         aln = nw.global_align(
#             self.query.sequence,
#             self.target.sequence,
#             gap_open=self._gap_open,
#             gap_extend=self._gap_extend,
#             matrix=matrix,
#         )
#         return aln

#     def _score_alignment(self):
#         if all([self._score_match is not None, self._score_mismatch is not None]):
#             matrix = self._get_matrix_file(
#                 match=self._score_match, mismatch=self._score_mismatch
#             )
#         elif self._matrix is not None:
#             matrix = self._get_matrix_file(matrix=self._matrix)
#         else:
#             matrix = self._get_matrix_file(match=self._match, mismatch=self._mismatch)
#         gap_open = (
#             self._score_gap_open if self._score_gap_open is not None else self._gap_open
#         )
#         gap_extend = (
#             self._score_gap_extend
#             if self._score_gap_extend is not None
#             else self._gap_extend
#         )
#         aln = nw.score_alignment(
#             self.aligned_query,
#             self.aligned_target,
#             gap_open=gap_open,
#             gap_extend=gap_extend,
#             matrix=matrix,
#         )
#         return aln

#     @staticmethod
#     def _build_matrix_from_params(match, mismatch, matrix_file):
#         mstring = " {}".format(match) if len(str(match)) == 1 else str(match)
#         mmstring = " {}".format(mismatch) if len(str(mismatch)) == 1 else str(mismatch)
#         residues = [
#             "A",
#             "C",
#             "D",
#             "E",
#             "F",
#             "G",
#             "H",
#             "I",
#             "K",
#             "L",
#             "M",
#             "N",
#             "P",
#             "Q",
#             "R",
#             "S",
#             "T",
#             "V",
#             "W",
#             "Y",
#             "*",
#         ]
#         header = "   " + "  ".join(residues)
#         matlist = [
#             header,
#         ]
#         for r1 in residues:
#             resline = [
#                 r1,
#             ]
#             for r2 in residues:
#                 resline.append(mstring if r1 == r2 else mmstring)
#             matlist.append(" ".join(resline))
#         open(matrix_file, "w").write("\n".join(matlist))
#         return matrix_file

#     @staticmethod
#     def _build_matrix_from_dict(matrix, matrix_file):
#         residues = sorted(matrix.keys())
#         header = "   " + "  ".join(residues)
#         matlist = [
#             header,
#         ]
#         for r1 in residues:
#             resline = [
#                 r1,
#             ]
#             for r2 in residues:
#                 s = matrix[r1][r2]
#                 score = " {}".format(s) if len(str(s)) == 1 else str(s)
#                 resline.append(score)
#             matlist.append(" ".join(resline))
#         with open(matrix_file, "w") as f:
#             f.write("\n".join(matlist))
#         return matrix_file

#     @staticmethod
#     def _get_builtin_matrix(matrix_name):
#         matrix_dir = os.path.join(
#             os.path.dirname(os.path.abspath(__file__)), "matrices"
#         )
#         matrices = [os.path.basename(f) for f in list_files(matrix_dir)]
#         if matrix_name.lower() not in matrices:
#             err = "The maxtrix name you provided ({}) is not built-in.".format(
#                 matrix_name
#             )
#             err += "Built in matrices are: {}".format(", ".join(matrices))
#             raise RuntimeError()
#         return os.path.join(matrix_dir, matrix_name.lower())


# -------------------------------------
#
#       VISUAL ALIGNMENT
#
# -------------------------------------


def dot_alignment(
    sequences,
    seq_field=None,
    name_field=None,
    root=None,
    root_name=None,
    cluster_threshold=0.75,
    as_fasta=False,
    just_alignment=False,
):
    """
    Creates a dot alignment (dots indicate identity, mismatches are represented by the mismatched
    residue) for a list of sequences.

    Args:

        sequence (list(Sequence)): A list of Sequence objects to be aligned.

        seq_field (str): Name of the sequence field key. Default is ``vdj_nt``.

        name_field (str): Name of the name field key. Default is ``seq_id``.

        root (str, Sequence): The sequence used to 'root' the alignment. This sequence will be at the
            top of the alignment and is the sequence against which dots (identity) will be evaluated.
            Can be provided either as a string corresponding to the name of one of the sequences in
            ``sequences`` or as a Sequence object. If not provided, ``sequences`` will be clustered
            at ``cluster_threshold`` and the centroid of the largest cluster will be used.

        root_name (str): Name of the root sequence. If not provided, the existing name of the root
            sequence (``name_field``) will be used. If ``root`` is not provided, the default ``root_name``
            is ``'centroid'``.

        cluster_threshold (float): Threshold with which to cluster sequences if ``root`` is not provided.
            Default is ``0.75``.

        as_fasta (bool): If ``True``, returns the dot alignment as a FASTA-formatted string, rather than
            a string formatted for human readability.

        just_alignment (bool): If ``True``, returns just the dot-aligned sequences as a list.

    Returns:

        If ``just_alignment`` is ``True``, a list of dot-aligned sequences (without sequence names) will be returned.
        If ``as_fasta`` is ``True``, a string containing the dot-aligned sequences in FASTA format will be returned.
        Otherwise, a formatted string containing the aligned sequences (with sequence names) will be returned.
    """
    import abstar

    from .cluster import cluster

    sequences = deepcopy(sequences)
    root = copy(root)

    # if custom seq_field is specified, copy to the .alignment_sequence attribute
    if seq_field is not None:
        if not all([seq_field in list(s.annotations.keys()) for s in sequences]):
            print(
                "\nERROR: {} is not present in all of the supplied sequences.\n".format(
                    seq_field
                )
            )
            sys.exit(1)
        for s in sequences:
            s.alignment_sequence = s[seq_field]
    else:
        for s in sequences:
            s.alignment_sequence = s.sequence

    # if custom name_field is specified, copy to the .id attribute
    if name_field is not None:
        if not all([name_field in list(s.annotations.keys()) for s in sequences]):
            print(
                "\nERROR: {} is not present in all of the supplied sequences.\n".format(
                    name_field
                )
            )
            sys.exit(1)
        for s in sequences:
            s.alignment_id = s[name_field]
    else:
        for s in sequences:
            s.alignment_id = s.id

    # parse the root sequence
    if all([root is None, root_name is None]):
        clusters = cluster(sequences, threshold=cluster_threshold, quiet=True)
        clusters = sorted(clusters, key=lambda x: x.size, reverse=True)
        centroid = clusters[0].centroid
        root = abstar.run(("centroid", centroid.sequence))
        root.alignment_id = "centroid"
        root.alignment_sequence = root[seq_field]
    elif type(root) in STR_TYPES:
        root = [s for s in sequences if s.alignment_id == root][0]
        if not root:
            print(
                "\nERROR: The name of the root sequence ({}) was not found in the list of input sequences.".format(
                    root
                )
            )
            print("\n")
            sys.exit(1)
        sequences = [s for s in sequences if s.alignment_id != root.alignment_id]
    elif type(root) == Sequence:
        if seq_field is not None:
            if seq_field not in list(root.annotations.keys()):
                print(
                    "\nERROR: {} is not present in the supplied root sequence.\n".format(
                        seq_field
                    )
                )
                sys.exit(1)
            root.alignment_sequence = root[seq_field]
        if name_field is not None:
            if name_field not in list(root.annotations.keys()):
                print(
                    "\nERROR: {} is not present in the supplied root sequence.\n".format(
                        name_field
                    )
                )
                sys.exit(1)
            root.alignment_id = root[name_field]
        sequences = [s for s in sequences if s.alignment_id != root.alignment_id]
    else:
        print(
            "\nERROR: If root is provided, it must be the name of a sequence \
              found in the supplied list of sequences or it must be a Sequence object."
        )
        print("\n")
        sys.exit(1)

    if root_name is not None:
        root.alignment_id = root_name
    else:
        root_name = root.alignment_id

    # compute and parse the alignment
    seqs = [(root.alignment_id, root.alignment_sequence)]
    seqs += [(s.alignment_id, s.alignment_sequence) for s in sequences]
    aln = muscle(seqs)
    g_aln = [a for a in aln if a.id == root_name][0]
    dots = [
        (root_name, str(g_aln.seq)),
    ]
    for seq in [a for a in aln if a.id != root_name]:
        s_aln = ""
        for g, q in zip(str(g_aln.seq), str(seq.seq)):
            if g == q == "-":
                s_aln += "-"
            elif g == q:
                s_aln += "."
            else:
                s_aln += q
        dots.append((seq.id, s_aln))
    if just_alignment:
        return [d[1] for d in dots]
    name_len = max([len(d[0]) for d in dots]) + 2
    dot_aln = []
    for d in dots:
        if as_fasta:
            dot_aln.append(">{}\n{}".format(d[0], d[1]))
        else:
            spaces = name_len - len(d[0])
            dot_aln.append(d[0] + " " * spaces + d[1])
    return "\n".join(dot_aln)
