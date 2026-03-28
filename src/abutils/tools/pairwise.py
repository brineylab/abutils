#!/usr/bin/env python
# filename: pairwise.py
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


import os
import sys
from abc import ABC
from collections.abc import Callable, Iterable
from itertools import groupby

import parasail
from Bio.SeqRecord import SeqRecord

from ..core.sequence import Sequence
from ..utils.decorators import lazy_property

__all__ = [
    "local_alignment",
    "global_alignment",
    "semiglobal_alignment",
    "PairwiseAlignment",
    "LocalAlignment",
    "GlobalAlignment",
    "SemiGlobalAlignment",
    "CIGAR",
]


# ----------------------------
#
#     PAIRWISE ALIGNMENT
#
# ----------------------------


class PairwiseAlignment(ABC):
    """
    Base class for local, global, and semi-global pairwise alignments.

    .. note::

        All comparisons between ``PairwiseAlignment`` objects are performed using the
        ``score`` attribute. This was done so that sorting alignments
        like so:

        .. code-block:: python

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
        matrix: parasail.Matrix | None = None,
        gap_open_penalty: int | None = None,
        gap_extend_penalty: int | None = None,
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
    def cigar(self) -> "CIGAR":
        """
        Alignment CIGAR, as a ``CIGAR`` object.

        Returns
        -------
        CIGAR
            CIGAR object representing the alignment.
        """
        cigar_string = self.alignment.cigar.decode.decode("utf-8")
        return CIGAR(cigar_string=cigar_string)

    @property
    def traceback(self):
        return self.alignment.traceback

    @lazy_property
    def score(self) -> int:
        """
        Alignment score.

        Returns
        -------
        int
            Alignment score.
        """
        return self.alignment.score

    @lazy_property
    def aligned_query(self) -> str:
        """
        Aligned query sequence.
        """
        return self.traceback.query

    @lazy_property
    def aligned_target(self) -> str:
        """
        Aligned target sequence.
        """
        return self.traceback.ref

    @lazy_property
    def alignment_midline(self) -> str:
        """
        Alignment midline. Pipe characters (|) indicate matches, and spaces indicate mismatches.
        """
        return self.traceback.comp.replace(".", " ")

    @lazy_property
    def query_begin(self) -> int:
        """
        Start position of the query sequence in the alignment.
        """
        query_begin = self.alignment.cigar.beg_query
        if self.cigar[0].element == "I":
            query_begin += self.cigar[0].length
        return query_begin

    @lazy_property
    def query_end(self) -> int:
        """
        End position of the query sequence in the alignment.
        """
        return self.alignment.end_query

    @lazy_property
    def target_begin(self) -> int:
        """
        Start position of the target sequence in the alignment.
        """
        target_begin = self.alignment.cigar.beg_ref
        if self.cigar[0].element == "D":
            target_begin += self.cigar[0].length
        return target_begin

    @lazy_property
    def target_end(self) -> int:
        """
        End position of the target sequence in the alignment.
        """
        return self.alignment.end_ref

    @property
    def target_id(self) -> str:
        """
        Sequence ID of the target sequence.
        """
        return self.target.id

    def align(self):
        # parasail throws an error if query or target sequences are of length 0
        # but -- annoyingly -- it also prints a warning to stdout
        # the stdout printing can mess with downstream code like abstar
        # to fix, we'll catch the problem before alignment and raise an error
        if len(self.query.sequence) == 0:
            raise ValueError(
                f"{self.alignment_type.upper()} ALIGNMENT ERROR: query sequence is empty."
            )
        if len(self.target.sequence) == 0:
            raise ValueError(
                f"{self.alignment_type.upper()} ALIGNMENT ERROR: target sequence is empty."
            )

        # do the alignment
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
        matrix: str | parasail.Matrix | None = None,
        gap_open_penalty: int | None = None,
        gap_extend_penalty: int | None = None,
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
        matrix: str | parasail.Matrix | None = None,
        gap_open_penalty: int | None = None,
        gap_extend_penalty: int | None = None,
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
        matrix: str | parasail.Matrix | None = None,
        gap_open_penalty: int | None = None,
        gap_extend_penalty: int | None = None,
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
        """
        Returns a list of CIGAR elements parsed from the CIGAR string. If the CIGAR string has not yet been parsed,
        it is parsed and the resulting list is cached for future use.

        :return: A list of `CIGARElement` objects representing the parsed CIGAR string.
        """
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
    query: str | SeqRecord | Sequence | Iterable,
    target: str | SeqRecord | Sequence | Iterable | None = None,
    targets: Iterable | None = None,
    match: int = 3,
    mismatch: int = -2,
    gap_open: int = -5,
    gap_extend: int = -2,
    matrix: str | parasail.Matrix | None = None,
    alignment_function: Callable = parasail.sw_trace_striped_16,
    aa: bool = False,
    gap_open_penalty: int | None = None,
    gap_extend_penalty: int | None = None,
) -> LocalAlignment | Iterable:
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
            - the name of a ``parasail`` `built-in sbustitution matrix <https://github.com/jeffdaily/parasail-python#substitution-matrices>`_.
            - path to a matrix file (in a format accepted by ``parasail``).
            - a ``parasail.Matrix`` object
        If not provided, a matrix will be created using `match` and `mismatch`.

    alignment_function : Callable, default=parasail.sw_trace_striped_16
        ``parasail`` `local alignment function <https://github.com/jeffdaily/parasail-python#standard-function-naming-convention>`_ to be used for alignment.

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
    query: str | SeqRecord | Sequence | Iterable,
    target: str | SeqRecord | Sequence | Iterable | None = None,
    targets: Iterable | None = None,
    match: int = 3,
    mismatch: int = -2,
    gap_open: int = -5,
    gap_extend: int = -2,
    matrix: str | parasail.Matrix | None = None,
    alignment_function: Callable = parasail.nw_trace_striped_16,
    aa: bool = False,
    gap_open_penalty: int | None = None,
    gap_extend_penalty: int | None = None,
) -> GlobalAlignment | Iterable:
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
            - the name of a ``parasail`` `built-in substitution matrix <https://github.com/jeffdaily/parasail-python#substitution-matrices>`_.
            - path to a matrix file (in a format accepted by ``parasail``).
            - a ``parasail.Matrix`` object
        If not provided, a matrix will be created using `match` and `mismatch`.

    alignment_function : Callable, default=parasail.nw_trace_striped_16
        ``parasail`` `global alignment function <https://github.com/jeffdaily/parasail-python#standard-function-naming-convention>`_ to be used for alignment.

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
    query: str | SeqRecord | Sequence | Iterable,
    target: str | SeqRecord | Sequence | Iterable | None = None,
    targets: Iterable | None = None,
    match: int = 3,
    mismatch: int = -2,
    gap_open: int = -5,
    gap_extend: int = -2,
    matrix: str | parasail.Matrix | None = None,
    alignment_function: Callable = parasail.sg_trace_striped_16,
    aa: bool = False,
    gap_open_penalty: int | None = None,
    gap_extend_penalty: int | None = None,
) -> SemiGlobalAlignment | Iterable:
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
            - the name of a ``parasail`` `built-in substitution matrix <https://github.com/jeffdaily/parasail-python#substitution-matrices>`_.
            - path to a matrix file (in a format accepted by ``parasail``).
            - a ``parasail.Matrix`` object
        If not provided, a matrix will be created using `match` and `mismatch`.

    alignment_function : Callable, default=parasail.sw_trace_striped_16
        ``parasail`` `semi-global alignment function <https://github.com/jeffdaily/parasail-python#standard-function-naming-convention>`_ to be used for alignment.

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
