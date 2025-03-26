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


import os
import platform
import random
import subprocess as sp
import sys
import tempfile
import uuid
from abc import ABC
from copy import copy, deepcopy
from io import StringIO
from itertools import groupby
from typing import Callable, Iterable, Optional, Union

import parasail
import pyfamsa
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..bin import get_path as get_binary_path
from ..core.sequence import Sequence
from ..io import read_fasta, to_fasta
from ..utils.decorators import lazy_property

__all__ = [
    "famsa",
    "mafft",
    "muscle",
    "muscle_v3",
    "local_alignment",
    "global_alignment",
    "semiglobal_alignment",
    "dot_alignment",
    "make_consensus",
    "MultipleSequenceAlignment",
    "PairwiseAlignment",
    "LocalAlignment",
    "GlobalAlignment",
    "SemiGlobalAlignment",
    "CIGAR",
]


# -------------------------------------
#
#     MULTIPLE SEQUENCE ALIGNMENT
#
# -------------------------------------


class MultipleSequenceAlignment:
    def __init__(
        self,
        input_alignment: Union[str, Iterable, MultipleSeqAlignment],
        fmt: str = "fasta",
    ):
        """
        Class for working with multiple sequence alignments.

        Parameters
        ----------
        input_alignment : (str, iterable, MultipleSeqAlignment)
            Can be one of several things:
                1. path to an alignment file
                2. an alignment string (for example, the result of calling ``read()`` on an alignment file)
                3. a BioPython ``MultipleSeqAlignment`` object
                4. a list of aligned BioPython ``SeqRecord`` objects
                5. a list of aligned ``abutils.Sequence`` objects

        fmt : str, default='fasta'
            Format of the input alignment. Choices are 'fasta', 'phylip', and 'clustal'.

        """
        self._sequences = None

        self.input = input_alignment
        self.fmt = fmt
        self.aln = self._process_alignment(input_alignment)

    def __getitem__(self, index):
        """Access part of the alignment. Drawn heavily from `BioPython's
        MultipleSeqAlignment object`_.

        Depending on the indices, you can get a Sequence object
        (representing a single row), a string (for a single character
        or a single column) or another MultipleSequenceAlignment object
        (representing some part or all of the alignment).

        - align[r,c] gives a single character as a string
        - align[r] gives a row as a Sequence
        - align[:,c] gives a column as a string

        align[:] and align[:,:] give a copy of the alignment

        Anything else gives a sub alignment, e.g.
        - align[0:2] or align[0:2,:] uses only row 0 and 1
        - align[:,1:3] uses only columns 1 and 2
        - align[0:2,1:3] uses only rows 0 & 1 and only cols 1 & 2

        We'll use the following example alignment here for illustration:

        .. code-block:: python
            from abutils import Sequence

            a = Sequence("AAAACGT", id="Alpha")
            b = Sequence("AAA-CGT", id="Beta")
            c = Sequence("AAAAGGT", id="Gamma")
            d = Sequence("AAAACGT", id="Delta")
            e = Sequence("AAA-GGT", id="Epsilon")

            align = MultipleSequenceAlignment([a, b, c, d, e])

        You can access a row of the alignment as a Sequence using an integer
        index (think of the alignment as a list of Sequence objects here):

        .. code-block:: python
            first_record = align[0]
            print("{first_record.id} {first_record.sequence}")
            # Alpha AAAACGT
            last_record = align[-1]
            print("{last_record.id} {last_record.sequence}")
            # Epsilon AAA-GGT

        You can also access use python's slice notation to create a sub-alignment
        containing only some of the Sequence objects:

        .. code-block:: python
            sub_alignment = align[2:5]
            print(sub_alignment)
            # Alignment with 3 rows and 7 columns
            # AAAAGGT Gamma
            # AAAACGT Delta
            # AAA-GGT Epsilon

        This includes support for a step, i.e. align[start:end:step], which
        can be used to select every second sequence:

        .. code-block:: python
            sub_alignment = align[::2]
            print(sub_alignment)
            # Alignment with 3 rows and 7 columns
            # AAAACGT Alpha
            # AAAAGGT Gamma
            # AAA-GGT Epsilon

        Or to get a copy of the alignment with the rows in reverse order:

        .. code-block:: python
            rev_alignment = align[::-1]
            print(rev_alignment)
            # Alignment with 5 rows and 7 columns
            # AAA-GGT Epsilon
            # AAAACGT Delta
            # AAAAGGT Gamma
            # AAA-CGT Beta
            # AAAACGT Alpha

        You can also use two indices to specify both rows and columns. Using simple
        integers gives you the entry as a single character string. e.g.

        .. code-block:: python
            align[3, 4]
            # 'C'

        This is equivalent to:

        .. code-block:: python
            align[3][4]
            # 'C'

        or:

        .. code-block:: python
            align[3].seq[4]
            # 'C'

        To get a single column (as a string) use this syntax:

        .. code-block:: python
            align[:, 4]
            # 'CCGCG'

        Or, to get part of a column,

        .. code-block:: python
            align[1:3, 4]
            # 'CG'

        However, in general you get a sub-alignment,

        .. code-block:: python
            print(align[1:5, 3:6])
            # Alignment with 4 rows and 3 columns
            # -CG Beta
            # AGG Gamma
            # ACG Delta
            # -GG Epsilon

        This should all seem familiar to anyone who has used the NumPy
        array or matrix objects.


        .. _BioPython's MultipleSeqAlignment object:
            https://biopython.org/DIST/docs/api/Bio.Align.MultipleSeqAlignment-class.html
        """
        # single indexing (e.g. align[6] or align[1:4])
        if isinstance(index, int):
            # e.g. result = align[x]
            # Return a single Sequence
            return self.sequences[index]
        elif isinstance(index, slice):
            # e.g. sub_align = align[i:j:k]
            return MultipleSequenceAlignment(self.sequences[index])
        elif len(index) != 2:
            raise TypeError("Invalid index type.")

        # double indexing (e.g. align[6, 1] or align[1:4, 6:10])
        row_index, col_index = index
        if isinstance(row_index, int):
            # e.g. row_or_part_row = align[6, 1:4], gives a Sequence
            return Sequence(
                self.sequences[row_index][col_index], id=self.sequences[row_index].id
            )
        elif isinstance(col_index, int):
            # e.g. col_or_part_col = align[1:5, 6], gives a string
            return "".join(seq[col_index] for seq in self.sequences[row_index])
        else:
            # e.g. sub_align = align[1:4, 5:7], gives a new alignment
            return MultipleSequenceAlignment(
                [
                    Sequence(seq[col_index], id=seq.id)
                    for seq in self.sequences[row_index]
                ]
            )

    def __iter__(self) -> Iterable[Sequence]:
        """
        Iterate over the aligned sequences in the alignment.
        """
        for s in self.sequences:
            yield s

    def __len__(self) -> int:
        """
        Return the number of sequences in the alignment.
        """
        return len(self.aln)

    def __repr__(self) -> str:
        """
        Return a string representation of the alignment.
        """
        return str(self.aln)

    def __str__(self) -> str:
        """
        Return a string representation of the alignment.
        """
        return self.aln_string

    def _process_alignment(self, input_alignment):
        if isinstance(input_alignment, MultipleSeqAlignment):
            return input_alignment
        elif isinstance(input_alignment, str):
            if os.path.isfile(input_alignment):
                return self._read_alignment_file(input_alignment, self.fmt)
            else:
                return self._read_alignment_string(input_alignment, self.fmt)
        elif isinstance(input_alignment, (list, tuple)):
            input_seqs = [Sequence(i) for i in input_alignment]
            aln_seqs = [SeqRecord(Seq(s.sequence), id=s.id) for s in input_seqs]
            aln = MultipleSeqAlignment(aln_seqs)
            self._sequences = input_seqs
            return aln
        else:
            raise ValueError(
                "Invalid input. Must be a path to an alignment file, a biopython \
                    ``MultipleSeqAlignment`` object, an alignment string, or a list of aligned Sequences."
            )

    def _read_alignment_file(self, aln_file, fmt):
        with open(aln_file, "r") as f:
            aln = AlignIO.read(f, fmt)
        return aln

    def _read_alignment_string(self, aln_file, fmt):
        return AlignIO.read(StringIO(aln_file), fmt)

    @property
    def aln_string(self) -> str:
        """
        Returns the alignment as a string.

        """
        return str(self.aln)

    @property
    def sequences(self) -> Iterable[Sequence]:
        """
        Sequences in the alignment. Note that this is a list of ``Sequence``
        objects, not BioPython ``SeqRecord`` objects, and that the sequences
        are aligned, meaning they are all the same length and may have gaps.

        Returns
        -------
        list of ``Sequence`` objects

        """
        if self._sequences is None:
            self._sequences = [Sequence(s) for s in self.aln]
        return self._sequences

    def format(self, format) -> str:
        """
        Format the alignment as a string.

        Parameters
        ----------
        format : str
            Format of the output alignment. Choices are 'fasta', 'phylip', and 'clustal'.

        """
        return self.aln.format(format)

    def make_consensus(
        self,
        name: Optional[str] = None,
        threshold: float = 0.51,
        ambiguous: Optional[str] = None,
        as_string: bool = False,
    ) -> Union[str, Sequence]:
        """
        Make a consensus sequence from the multiple sequence alignment.

        Parameters
        ----------
        name : str, optional
            Name of the consensus sequence. If not provided, a random UUID
            will be used.

        threshold : float, default=0.51
            Threshold for calling a consensus base. Default is 0.51, meaning
            that a base will be called if it is present in at least 51% of
            the sequences.

        ambiguous : str, optional
            Character to use for ambiguous bases. If not provided, the
            default is "N" if the sequences contain only standard nucleotides
            (A, T, G, C), and "X" otherwise.

        as_string : bool, default=False
            If ``True``, returns the consensus sequence as a string. If
            ``False``, returns a ``Sequence`` object.

        Returns
        -------
        consensus : str or ``Sequence``
            If ``as_string`` is ``True``, returns the consensus sequence as
            a string. If ``False``, returns a ``Sequence`` object.

        """
        if ambiguous is None:
            if all(
                [
                    nt in ["A", "T", "G", "C", "N", "-"]
                    for s in self.sequences
                    for nt in s.sequence
                ]
            ):
                ambiguous = "N"
            else:
                ambiguous = "X"
        summary_align = AlignInfo.SummaryInfo(self.aln)
        consensus = summary_align.gap_consensus(
            threshold=threshold, ambiguous=ambiguous
        )
        consensus_string = str(consensus).replace("-", "")
        if as_string:
            return consensus_string
        return Sequence(
            consensus_string, id=name if name is not None else str(uuid.uuid4())
        )


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
        Path to a FASTA-formatted input file. *Depricated* (use `sequences`
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
    ffile = to_fasta(
        sequences, id_key=id_key, sequence_key=seq_key, tempfile_dir="/tmp"
    )
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
        mafft_bin = get_binary_path("mafft")
        # mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # # mafft_bin = os.path.join(
        # #     mod_dir, "bin/mafft_{}".format(platform.system().lower())
        # # )
        # system = platform.system().lower()
        # machine = platform.machine().lower().replace("x86_64", "amd64")
        # mafft_bin = os.path.join(mod_dir, f"bin/mafft_{system}_{machine}/mafft.bat")
        # # if system == "darwin":
        # #     mafft_bin = os.path.join(mod_dir, f"bin/mafft_{system}_{machine}/mafft.bat")
        # # else:
        # #     mafft_bin = os.path.join(mod_dir, f"bin/mafft_{system}_{machine}")
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
        print("mafft binary path:", mafft_bin)
        print("mafft command:", mafft_cline)
        print(stdout)
        print(stderr)
    if os.stat(alignment_file).st_size == 0:
        err = "WARNING: output alignment file is empty. "
        err += "Verify that MAFFT is installed and the input data is valid,"
        print(err)
        return None
    if as_file:
        return alignment_file

    if as_string:
        with open(alignment_file, "r") as f:
            aln = f.read()
    else:
        # aln = AlignIO.read(f, fmt)
        aln = MultipleSequenceAlignment(alignment_file, fmt.lower())
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
    # >>>>>> TEMPORARY FIX for MUSCLE 5.0 and MacOS
    if platform.system().lower() == "darwin" and muscle_bin is None:
        return muscle_v3(
            sequences=sequences,
            alignment_file=alignment_file,
            as_file=as_file,
            as_string=as_string,
            fmt="fasta",
            id_key=id_key,
            seq_key=seq_key,
            debug=debug,
            fasta=fasta,
        )
    # <<<<<< END TEMPORARY FIX

    # process input
    if fasta is not None:
        sequences = fasta
    ffile = to_fasta(
        sequences, id_key=id_key, sequence_key=seq_key, tempfile_dir="/tmp"
    )
    # configure output path
    if alignment_file is None:
        # as_file = False
        alignment_file = tempfile.NamedTemporaryFile(delete=False).name
    else:
        alignment_file = os.path.abspath(alignment_file)
    # muscle binary
    if muscle_bin is None:
        muscle_bin = get_binary_path("muscle")
        # mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # system = platform.system().lower()
        # machine = platform.machine().lower().replace("x86_64", "amd64")
        # muscle_bin = os.path.join(mod_dir, f"bin/muscle_{system}_{machine}")
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
        err += "Verify that MUSCLE is installed and the input data is valid"
        print(err)
        return None
    if as_file:
        return alignment_file
    if as_string:
        with open(alignment_file, "r") as f:
            aln = f.read()
    else:
        # aln = AlignIO.read(f, "fasta")
        aln = MultipleSequenceAlignment(alignment_file, "fasta")
    os.unlink(alignment_file)
    return aln


def famsa(
    sequences: Union[str, Iterable],
    alignment_file: Optional[str] = None,
    fmt: str = "fasta",
    as_file: bool = False,
    as_string: bool = False,
    threads: int = 0,
    guide_tree: str = "sl",
    tree_heuristic: Optional[str] = None,
    medoid_threshold: int = 0,
    n_refinements: int = 100,
    keep_duplicates: bool = False,
    refine: Optional[bool] = None,
    id_key: Optional[str] = None,
    seq_key: Optional[str] = None,
    debug: bool = False,
    fasta: Optional[str] = None,
) -> Union[MultipleSeqAlignment, str]:
    """
    Performs multiple sequence alignment with `FAMSA`_ using the `pyfamsa`_ package.

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

    threads : int, default=0
        Number of threads for FAMSA to use. Default is ``0``, which uses all
        available CPUs.

    guide_tree : str, default='sl'
        Method for building the guide tree. Choices are 'sl', 'slink', 'nj', and 'upgma'.

    tree_heuristic : str, default=None
        The heuristic to use for constructing the tree. Supported values are:
        'medoid', 'part', or ``None`` to disable heuristics.

    medoid_threshold : int, default=0
        Minimum number of sequences a set must contain for medoid trees to be used,
        if enabled with ``tree_heuristic``. Default is ``0``.

    n_refinements : int, default=100
        Number of refinement iterations to perform. Default is ``100``.

    keep_duplicates : bool, default=False
        If ``True``, duplicate sequences are kept in the alignment. Default is ``False``.

    refine : bool, default=None
        If ``True``, the alignment is refined. If ``False``, the alignment is not refined.
        If ``None``, the alignment is refined if the number of sequences is less than 1000.

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


    .. _FAMSA:
        https://github.com/refresh-bio/FAMSA
    .. _pyfamsa:
        https://github.com/althonos/pyfamsa
    """
    # process input
    if fasta is not None:
        sequences = fasta
    if isinstance(sequences, str):
        seqs = read_fasta(sequences)
    else:
        seqs = [Sequence(s, id_key=id_key, seq_key=seq_key) for s in sequences]
        # for s in sequences:
        #     _seq = s.sequence if seq_key is None else s[seq_key]
        #     _id = s.id if id_key is None else s[id_key]
        #     seqs.append(Sequence(_seq, id=_id))
    pyfamsa_seqs = [pyfamsa.Sequence(s.id.encode(), s.sequence.encode()) for s in seqs]
    # do the alignment
    aligner = pyfamsa.Aligner(
        threads=threads,
        guide_tree=guide_tree,
        tree_heuristic=tree_heuristic,
        medoid_threshold=medoid_threshold,
        n_refinements=n_refinements,
        keep_duplicates=keep_duplicates,
        refine=refine,
    )
    aln = aligner.align(pyfamsa_seqs)
    aligned_seqs = [Sequence(s.sequence.decode(), id=s.id.decode()) for s in aln]
    msa = MultipleSequenceAlignment(aligned_seqs)
    # build output
    if as_file:
        if alignment_file is None:
            alignment_file = tempfile.NamedTemporaryFile(delete=False).name
        with open(alignment_file, "w") as f:
            f.write(msa.format(fmt))
        return alignment_file
    elif as_string:
        return msa.format(fmt)
    return msa


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
    # if fasta is not None:
    #     sequences = fasta
    # fasta_string = to_fasta(
    #     sequences,
    #     as_string=True,
    #     id_key=id_key,
    #     sequence_key=seq_key,
    # )
    if fasta is not None:
        sequences = fasta
    ffile = to_fasta(
        sequences, id_key=id_key, sequence_key=seq_key, tempfile_dir="/tmp"
    )
    # configure output path
    if alignment_file is None:
        # as_file = False
        alignment_file = tempfile.NamedTemporaryFile(delete=False).name
    else:
        alignment_file = os.path.abspath(alignment_file)
    # configure output path
    if alignment_file is None:
        as_file = False
        alignment_file = tempfile.NamedTemporaryFile(delete=False).name
    else:
        alignment_file = os.path.abspath(alignment_file)

    # muscle binary
    if muscle_bin is None:
        muscle_bin = get_binary_path("muscle3")
        # mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # system = platform.system().lower()
        # muscle_bin = os.path.join(mod_dir, f"bin/muscle3_{system}")
    # do the alignment
    aln_lookup = {
        "fasta": "-fastaout",
        "phylip": "-phyiout",
        "clustal": "-clwstrict -clwout",
    }
    if fmt.lower() not in aln_lookup.keys():
        err = f"\tERROR: invalid alignment format: '{fmt.lower()}'\n"
        err += "Valid formats are 'fasta', 'phylip', and 'clustal'.\n"
        raise ValueError(err)
    aln_format = aln_lookup[fmt.lower()]
    muscle_cline = f"{muscle_bin} -in {ffile} {aln_format} {alignment_file}"
    if maxiters is not None:
        muscle_cline += f" -maxiters {maxiters}"
    if diags:
        muscle_cline += " -diags"
    if all([gap_open is not None, gap_extend is not None]):
        muscle_cline += f" -gapopen {gap_open} -gapextend {gap_extend}"
    if debug:
        print("muscle binary path:", muscle_bin)
        print("muscle command:", muscle_cline)
    muscle = sp.Popen(
        str(muscle_cline),
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        universal_newlines=True,
        shell=True,
    )
    stdout, stderr = muscle.communicate()
    if debug:
        print(stdout)
        print(stderr)
    # output
    if os.stat(alignment_file).st_size == 0:
        err = "WARNING: output alignment file is empty. "
        err += "Verify that MUSCLE is installed and the input data is valid"
        print(err)
        return None
    if as_file:
        return alignment_file
    if as_string:
        with open(alignment_file, "r") as f:
            aln = f.read()
    else:
        # aln = AlignIO.read(f, fmt)
        aln = MultipleSequenceAlignment(alignment_file, fmt.lower())
    os.unlink(alignment_file)
    return aln

    # if as_string:
    #     return alignment
    # elif as_file:
    #     with open(alignment_file, "w") as f:
    #         f.write(alignment)
    #     return alignment_file
    # else:
    #     return AlignIO.read(StringIO(alignment), fmt)


# def consensus(aln, name=None, threshold=0.51, ambiguous="N"):
#     summary_align = AlignInfo.SummaryInfo(aln)
#     consensus = summary_align.gap_consensus(threshold=threshold, ambiguous=ambiguous)
#     if name is None:
#         name = str(uuid.uuid4())
#     consensus_string = str(consensus).replace("-", "")
#     return (name, consensus_string.upper())


def make_consensus(
    sequences: Optional[Iterable[Sequence]] = None,
    alignment: Union[str, MultipleSeqAlignment, MultipleSequenceAlignment, None] = None,
    algo: str = "mafft",
    alignment_kwargs: Optional[dict] = None,
    name: Optional[str] = None,
    downsample_to: Optional[int] = None,
    threshold: float = 0.51,
    ambiguous: Optional[str] = None,
    seed: Optional[int] = None,
    as_string: bool = False,
) -> Sequence:
    """
    Makes a consensus sequence from a list of sequences.

    Parameters
    ----------
    sequences : iterable, optional
        List of sequences from which to make a consensus. Sequences should not already
        be aligned (if they are, use ``alignment`` instead). One of ``sequences``
        or ``alignment`` must be provided.

    alignment : str or ``MultipleSeqAlignment``, optional
        A FASTA-formatted alignment string, a path to a FASTA-formatted alignment file,
        a BioPython ``MultipleSeqAlignment`` object, or a ``MultipleSequenceAlignment`` object.
        One of ``sequences`` or ``alignment`` must be provided.

    algo : str, default='mafft'
        Algorithm to use for multiple sequence alignment. Choices are 'mafft',
        'famsa', or 'muscle'.

    alignment_kwargs : dict, optional
        Additional keyword arguments to pass to the alignment function.

    name : str, optional
        Name for the consensus sequence. If not provided, a random UUID will be used.

    downsample_to : int, optional
        If provided, downsamples the input sequences to the specified number.

    threshold : float, default=0.51
        Threshold for consensus sequence generation. Default is 0.51.

    ambiguous : str, default='N'

    Returns
    -------
    consensus : ``Sequence``
        Consensus sequence.

    """
    if sequences is None and alignment is None:
        raise ValueError("Must provide either <sequences> or <alignment>.")
    elif sequences is not None:
        # downsampling
        if downsample_to is not None and len(sequences) > downsample_to:
            if seed is not None:
                random.seed(seed)
            sequences = random.sample(sequences, downsample_to)
        # alignment
        if alignment_kwargs is None:
            alignment_kwargs = {}
        alignment_kwargs["as_file"] = False
        if algo.lower() == "mafft":
            aln = mafft(sequences, **alignment_kwargs)
        elif algo.lower() == "muscle":
            aln = muscle(sequences, **alignment_kwargs)
        elif algo.lower() == "famsa":
            aln = famsa(sequences, **alignment_kwargs)
        else:
            err = f"Invalid algorithm: {algo}. Must be 'mafft', 'famsa', or 'muscle'."
            raise ValueError(err)
    else:
        if isinstance(alignment, MultipleSequenceAlignment):
            aln = alignment
        else:
            aln = MultipleSequenceAlignment(alignment)
    # consensus
    return aln.make_consensus(
        name=name, threshold=threshold, ambiguous=ambiguous, as_string=as_string
    )
    # summary_align = AlignInfo.SummaryInfo(aln)
    # consensus = summary_align.gap_consensus(threshold=threshold, ambiguous=ambiguous)
    # consensus_string = str(consensus).replace("-", "").upper()
    # if as_string:
    #     return consensus_string
    # return Sequence(
    #     consensus_string, id=name if name is not None else str(uuid.uuid4())
    # )


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
        # return self.alignment.cigar.beg_query
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
        # return self.alignment.cigar.beg_ref
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
    elif isinstance(root, str):
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
