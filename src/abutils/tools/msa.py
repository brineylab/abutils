#!/usr/bin/env python
# filename: msa.py
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
from collections.abc import Iterable
from copy import copy, deepcopy
from io import StringIO

import pyfamsa
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..bin import get_path as get_binary_path
from ..core.sequence import Sequence
from ..io import read_fasta, to_fasta

__all__ = [
    "famsa",
    "mafft",
    "muscle",
    "muscle_v3",
    "make_consensus",
    "dot_alignment",
    "MultipleSequenceAlignment",
]


# -------------------------------------
#
#     MULTIPLE SEQUENCE ALIGNMENT
#
# -------------------------------------


class MultipleSequenceAlignment:
    def __init__(
        self,
        input_alignment: str | Iterable | MultipleSeqAlignment,
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
        name: str | None = None,
        threshold: float = 0.51,
        ambiguous: str | None = None,
        as_string: bool = False,
    ) -> str | Sequence:
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
        # Build consensus using position-wise frequency analysis
        # This replaces the deprecated gap_consensus method from Biopython
        consensus_chars = []
        alignment_length = self.aln.get_alignment_length()
        num_sequences = len(self.aln)
        for i in range(alignment_length):
            column = self.aln[:, i]
            # Count occurrences of each character
            char_counts = {}
            for char in column:
                char_counts[char] = char_counts.get(char, 0) + 1
            # Find the most common character meeting threshold
            best_char = ambiguous
            best_freq = 0
            for char, count in char_counts.items():
                freq = count / num_sequences
                if freq >= threshold and freq > best_freq:
                    best_char = char
                    best_freq = freq
            consensus_chars.append(best_char)
        consensus_string = "".join(consensus_chars).replace("-", "")
        if as_string:
            return consensus_string
        return Sequence(
            consensus_string, id=name if name is not None else str(uuid.uuid4())
        )


def mafft(
    sequences: str | Iterable,
    alignment_file: str | None = None,
    fmt: str = "fasta",
    threads: int = -1,
    as_file: bool = False,
    as_string: bool = False,
    reorder: bool = True,
    mafft_bin: str | None = None,
    id_key: str | None = None,
    seq_key: str | None = None,
    debug: bool = False,
    fasta: str | None = None,
) -> MultipleSeqAlignment | str:
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
        aln = MultipleSequenceAlignment(alignment_file, fmt.lower())
    os.unlink(alignment_file)
    return aln


def muscle(
    sequences: str | Iterable,
    alignment_file: str | None = None,
    as_file: bool = False,
    as_string: bool = False,
    muscle_bin: str | None = None,
    threads: int | None = None,
    id_key: str | None = None,
    seq_key: str | None = None,
    debug: bool = False,
    fasta: str | None = None,
) -> MultipleSeqAlignment | str:
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
        aln = MultipleSequenceAlignment(alignment_file, "fasta")
    os.unlink(alignment_file)
    return aln


def famsa(
    sequences: str | Iterable,
    alignment_file: str | None = None,
    fmt: str = "fasta",
    as_file: bool = False,
    as_string: bool = False,
    threads: int = 0,
    guide_tree: str = "sl",
    tree_heuristic: str | None = None,
    medoid_threshold: int = 0,
    n_refinements: int = 100,
    keep_duplicates: bool = False,
    refine: bool | None = None,
    id_key: str | None = None,
    seq_key: str | None = None,
    debug: bool = False,
    fasta: str | None = None,
) -> MultipleSeqAlignment | str:
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
    sequences: str | Iterable | None = None,
    alignment_file: str | None = None,
    fmt: str = "fasta",
    as_file: bool = False,
    as_string: bool = False,
    muscle_bin: str | None = None,
    maxiters: int | None = None,
    diags: int | None = None,
    gap_open: int | None = None,
    gap_extend: int | None = None,
    id_key: str | None = None,
    seq_key: str | None = None,
    debug: bool = False,
    fasta: str | None = None,
) -> MultipleSeqAlignment | str:
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
    # do the alignment
    aln_lookup = {
        "fasta": "-fastaout",
        "phylip": "-phyiout",
        "clustal": "-clwstrict -clwout",
    }
    if fmt.lower() not in aln_lookup:
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
        aln = MultipleSequenceAlignment(alignment_file, fmt.lower())
    os.unlink(alignment_file)
    return aln


def make_consensus(
    sequences: Iterable[Sequence] | None = None,
    alignment: str | MultipleSeqAlignment | MultipleSequenceAlignment | None = None,
    algo: str = "mafft",
    alignment_kwargs: dict | None = None,
    name: str | None = None,
    downsample_to: int | None = None,
    threshold: float = 0.51,
    ambiguous: str | None = None,
    seed: int | None = None,
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
