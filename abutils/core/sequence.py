#!/usr/bin/env python
# filename: sequence.py


#
# Copyright (c) 2015 Bryan Briney
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


import csv
import gzip
import json
import operator
import os
import sys
import tempfile
import uuid
from collections import OrderedDict
from typing import Iterable, Optional, Union

import dnachisel as dc
import pandas as pd
import polars as pl
import pyfastx
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..utils.codons import codon_lookup
from ..utils.path import make_dir
from ..utils.utilities import nested_dict_lookup

STR_TYPES = [str]


class Sequence(object):
    """
    Container for biological (RNA, DNA, or protein) sequences.

    ``seq`` can be one of several things:

        1) a raw sequence, as a string

        2) an iterable, formatted as ``[seq_id, sequence]``

        3) a dict, containing at least the sequence ID and a
           sequence. Alternate ``id_key`` and ``seq_key``
           can be provided at instantiation.

        4) a Biopython ``SeqRecord`` object

        5) an abutils ``Sequence`` object

    If ``seq`` is provided as a string, the sequence ID can optionally be
    provided via the ``id`` keyword argument.  If ``seq`` is a string and ``id`` is not provided,
    a random sequence ID will be generated with ``uuid.uuid4()``.

    Quality scores can be supplied with ``qual`` or as part of a ``SeqRecord`` object.
    If providing both a SeqRecord object with quality scores and quality scores
    via ``qual``, the ``qual`` scores will override the SeqRecord quality scores.

    If ``seq`` is a dictionary, typically the result of a MongoDB query, the dictionary
    can be accessed directly from the ``Sequence`` instance (via the ``annotations`` property). To retrive the value
    for ``'junc_aa'`` in the instantiating dictionary, you would simply:

    .. code-block:: python

        s = Sequence(dict)
        junc = s['junc_aa']

    If ``seq`` is a dictionary, an optional ``id_key`` and ``seq_key`` can be provided,
    which tells the ``Sequence`` object which field to use to populate ``Sequence.id`` and
    ``Sequence.sequence``. Defaults for both ``id_key`` and ``seq_key`` are ``None``, which
    results in abutils trying to determine the appropriate key. For ``id_key``, the following
    keys are tried: ``['seq_id', 'sequence_id']``. For ``seq_key``, the following keys are tried:
    ``['vdj_nt', 'sequence_nt', 'sequence']``. If none of the attempts are successful, the
    ``Sequence.id`` or ``Sequence.sequence`` attributes will be ``None``.

    Alternately, the ``__getitem__()`` interface can be used to obtain a slice from the
    ``sequence`` attribute. An example of the distinction:

    .. code-block:: python

        d = {'name': 'MySequence', 'sequence': 'ATGC'}
        seq = Sequence(d, id_key='name', seq_key='sequence')

        seq['name']  # 'MySequence'
        seq[:2]  # 'AT'

    If the ``Sequence`` is instantiated with a dictionary, calls to ``__contains__()`` will
    return ``True`` if the supplied item is a key in the dictionary. In non-dict instantiations,
    ``__contains__()`` will look in the ``Sequence.sequence`` field directly (essentially a
    motif search). For example:

    .. code-block:: python

        dict_seq = Sequence({'seq_id': 'seq1', 'vdj_nt': 'ACGT'})
        'seq_id' in dict_seq  # TRUE
        'ACG' in dict_seq     # FALSE

        str_seq = Sequence('ACGT', id='seq1')
        'seq_id' in str_seq  # FALSE
        'ACG' in str_seq     # TRUE

    .. note::

        When comparing ``Sequence`` objects, they are comsidered equal only if their
        sequences and IDs are identical. This means that two ``Sequence`` objects
        with identical sequences but without user-supplied IDs won't be equal,
        because their IDs will have been randomly generated.
    """

    def __init__(
        self,
        sequence: Union[str, Iterable, dict, SeqRecord],
        id: Optional[str] = None,
        qual: Optional[str] = None,
        annotations: Optional[dict] = None,
        id_key: Optional[str] = None,
        seq_key: Optional[str] = None,
    ):
        super().__init__()
        self._input_sequence = None
        self._input_id = id
        self._input_qual = qual
        self._annotations = annotations
        self._clobbered_annotations = {}
        self.id = None
        self.sequence = None
        self.qual = None
        self._fasta = None
        self._fastq = None
        self._reverse_complement = None
        self._codon_optimized = None
        self._strand = None
        self.id_key = id_key
        self.seq_key = seq_key
        self.description = None

        self._process_input(sequence, id, qual)

    def __len__(self) -> int:
        """
        int: Returns the length of the ``sequence`` attribute
        """
        return len(self.sequence)

    def __iter__(self) -> Iterable[str]:
        """iter: Returns an iterator over the ``sequence`` attribute"""
        return iter(self.sequence)

    def __reversed__(self) -> str:
        """
        str: Returns the reverse of the ``sequence`` attribute
        """
        return "".join(reversed(self.sequence))

    def __contains__(self, item: str) -> bool:
        """
        Returns True if the supplied item is in the ``Sequence`` object.

        Parameters
        ----------

        item : str
            Item to be checked for membership in the ``Sequence`` object.

        Returns
        -------
        bool
            If the instance was initialzed with a dictonary (which means
            the ``annotations`` attribute is not empty), ``__contains__(key)``
            will return ``True`` if ``key`` is in ``annotations.keys()``. If ``annotations``
            is an empty dict, indicating instantiation without a dictionary,
            ``__contains__(motif)`` will return True if ``motif`` is in the
            ``sequence`` attribute.

        """
        if self.annotations:
            return item in self.annotations.keys()
        return item in self.sequence

    def __getitem__(self, key) -> Union[str, dict, None]:
        """
        Returns a slice of the ``Sequence`` object.

        Parameters
        ----------
        key : Union[int, str, slice]
            If ``Sequence`` was instantiated with a dictionary, ``key`` can be a string
            that corresponds to a key in the dictionary. If ``Sequence`` was instantiated
            with a string, ``key`` can be an integer or a slice.

        Returns
        -------
        Union[str, dict, None]
            If ``key`` is a string and the ``Sequence`` was instantiated with a dictionary,
            the value for the key will be returned. If ``key`` is a string and the ``Sequence``
            was instantiated with a string, ``None`` will be returned. If ``key`` is an integer
            or a slice, the corresponding slice of the ``sequence`` attribute will be returned.

        """
        if type(key) == slice:
            return self.sequence[key]
        elif key in self.annotations.keys():
            return self.annotations.get(key, None)
        elif isinstance(key, int):
            return self.sequence[key]
        return None

    def __setitem__(
        self, key: Union[str, int, float], val: Union[str, int, float]
    ) -> None:
        """
        Sets a value in the ``annotations`` property of a ``Sequence`` object.

        Parameters
        ----------
        key : Union[str, int, float]
            Key for the annotation to be set.

        val : Union[str, int, float]
            Value to be set for the annotation.

        """
        self.annotations[key] = val

    def __eq__(self, other: "Sequence") -> bool:
        """
        Returns True if the sequences and IDs of two ``Sequence`` objects are identical.

        Parameters
        ----------
        other : Sequence
            Another ``Sequence`` object.

        """
        if all([hasattr(other, "sequence"), hasattr(other, "id")]):
            return all([self.sequence == other.sequence, self.id == other.id])
        return False

    # def __getattr__(self, name):
    #     """
    #     Allows attribute access to the ``annotations`` dictionary.

    #     Parameters
    #     ----------
    #     name : str
    #         Name of the attribute to be accessed.

    #     Returns
    #     -------
    #     Any
    #         The value of the attribute, if present in the ``annotations`` dictionary.

    #     Raises
    #     ------
    #     AttributeError
    #         If the attribute is not in the ``annotations`` dictionary.
    #     """
    #     if name in self.annotations:
    #         return self.annotations[name]
    #     raise AttributeError(f"{self.__class__.__name__} has no attribute '{name}'")

    # # def __setattr__(self, name, value):
    # #     """
    # #     Allows attribute assignment to the ``annotations`` dictionary.

    # #     Parameters
    # #     ----------
    # #     name : str
    # #         Name of the attribute to be set.

    # #     value : Any
    # #         Value to be set for the attribute.
    # #     """
    # #     # start with the built-in attributes
    # #     if name in self.__dict__:
    # #         super().__setattr__(name, value)
    # #     # then try the annotations
    # #     elif name in self._annotations:
    # #         self.annotations[name] = value
    # #     # otherwise, set the attribute normally (make a new attribute)
    # #     else:
    # #         super().__setattr__(name, value)

    @classmethod
    def from_json(cls, data):
        if isinstance(data, str):
            data = json.loads(data)
        s = cls(None)
        s.__dict__ = data
        return s

    @property
    def fasta(self) -> str:
        """
        Returns the sequence, as a FASTA-formatted string.

        .. note::

                The FASTA string is built using ``Sequence.id`` and ``Sequence.sequence``.

        Returns
        -------
        str
            Returns the sequence, as a FASTA-formatted string.

        """
        if not self._fasta:
            self._fasta = ">{}\n{}".format(self.id, self.sequence)
        return self._fasta

    @property
    def fastq(self):
        """
        Returns the sequence, as a FASTQ-formatted string.

        .. note::

            The FASTQ string is built using ``Sequence.id``, ``Sequence.sequence``, and
            ``Sequence.qual``.

        Returns
        -------
        str
            Returns the sequence, as a FASTQ-formatted string

        """
        if self.qual is None:
            self._fastq = None
        else:
            if self._fastq is None:
                self._fastq = "@{}\n{}\n+\n{}".format(self.id, self.sequence, self.qual)
        return self._fastq

    @property
    def reverse_complement(self) -> Union[str, "Sequence"]:
        """
        Returns the reverse complement of ``Sequence.sequence``.

        Returns
        -------
            Returns a ``str`` if ``in_place`` is ``False``, otherwise returns an updated
            ``Sequence`` object in which the ``sequence`` property has been replaced
             with the reverse complement.

        """
        return reverse_complement(self, in_place=False)
        # if self._reverse_complement is None:
        #     self._reverse_complement = self._get_reverse_complement()
        # return self._reverse_complement

    @property
    def codon_optimized(self):
        if self._codon_optimized is None:
            self._codon_optimized = self.codon_optimize(as_string=True)
        return self._codon_optimized

    @codon_optimized.setter
    def codon_optimized(self, val):
        self._codon_optimized = val

    @property
    def annotations(self):
        """
        Annotations is a dictionary that contains any additional information
        about the sequence. This can include sequence annotations, such as
        VDJ annotations, or any other information that might be useful to
        store with the sequence (like donor, group, etc).

        """
        if self._annotations is None:
            self._annotations = {}
        return self._annotations

    @annotations.setter
    def annotations(self, val):
        self._annotations = val

    @property
    def strand(self):
        if self._strand is None:
            self._strand = "plus"
        return self._strand

    @strand.setter
    def strand(self, strand):
        self._strand = strand

    def translate(self, sequence_key: Optional[str] = None, frame: int = 1) -> str:
        """
        Translate a nucleotide sequence.

        Parameters
        ----------
        sequence_key : str, default=None
            Name of the annotation field containg the sequence to be translated.
            If not provided, ``Sequence.sequence`` is used.

        frame : int, default=1
            Reading frame to translate. Default is ``1``.

        Returns
        -------
        str
            Translated sequence

        """
        return translate(self, sequence_key, frame)

    def codon_optimize(
        self,
        sequence_key: str = "sequence_aa",
        id_key: str = "sequence_id",
        frame: Optional[int] = None,
        as_string: bool = True,
    ) -> Union[str, "Sequence"]:
        """
        Codon optimize a sequence.

        Parameters
        ----------
        sequence_key : str, default="sequence_aa"
            Name of the annotation field containg the sequence to be translated.
            If not provided, ``Sequence.sequence`` is used.

        id_key : str, default="sequence_id"
            Name of the annotation field containg the sequence id.
            If not provided, ``Sequence.id`` is used.

        frame : int, default=1
            Reading frame to translate. Default is ``1``.

        as_string : bool, default=True
            If ``True``, the optimized sequence will be returned as a ``str``.
            If ``False``, the optimized sequence will be returned as a ``Sequence`` object.

        Returns
        -------
        Union[str, Sequence]
            Translated sequence as a ``str`` (if ``as_string`` is ``True``) or
            ``Sequence`` object (if ``as_string`` is ``False``).

        """
        return codon_optimize(
            sequence=self,
            sequence_key=sequence_key,
            id_key=id_key,
            frame=frame,
            as_string=as_string,
        )

    def as_fasta(
        self, name_field: Optional[str] = None, seq_field: Optional[str] = None
    ) -> str:
        """
        Returns the sequence, as a FASTA-formatted string.

        Parameters
        ----------
        name_field : str, default=None
            Name of the annotation field containing the sequence name. If not provided,
            ``Sequence.id`` is used.

        seq_field : str, default=None
            Name of the annotation field containing the sequence. If not provided,
            ``Sequence.sequence`` is used.

        Returns
        -------
        str
            Returns the sequence, as a FASTA-formatted string.

        """
        name = None
        sequence = None
        if name_field is not None:
            name = self.annotations[name_field]
        if seq_field is not None:
            sequence = self.annotations[seq_field]
        if name is None:
            name = self.id
        if sequence is None:
            sequence = self.sequence
        return ">{}\n{}".format(name, sequence)

    def region(self, start=0, end=None):
        """
        Returns a region of ``Sequence.sequence``, in FASTA format.

        If called without kwargs, the entire sequence will be returned.

        Parameters
        ----------
        start : int, default=0
            Start position of the region to be returned. Default is 0.

        end : int, default=None
            End position of the region to be returned. Negative values

        Returns
        -------
        str
            A region of ``Sequence.sequence``, in FASTA format

        """
        if end is None:
            end = len(self.sequence)
        return ">{}\n{}".format(self.id, self.sequence[start:end])

    def keys(self):
        """
        Returns the keys of the ``annotations`` attribute.
        """
        return self.annotations.keys()

    def values(self):
        """
        Returns the values of the ``annotations`` attribute.
        """
        return self.annotations.values()

    def get(
        self, key: Union[str, int, float], default: Union[str, int, float, None] = None
    ) -> Union[str, int, float, None]:
        """
        Returns the value of a key in the ``annotations`` attribute.

        Parameters
        ----------
        key : Union[str, int, float]
            Key for the annotation to be returned.

        default : Union[str, int, float, None], default=None
            Value to be returned if the key is not in the ``annotations`` attribute.

        Returns
        -------
        Union[str, int, float, None]
            Value of the key in the ``annotations`` attribute.

        """
        return self.annotations.get(key, default)

    def _process_input(self, seq, id, qual):
        # sequence was passed as a string
        if type(seq) in STR_TYPES:
            if id is None:
                id = str(uuid.uuid4())
            self.sequence = str(seq).upper()
            self.id = id
            self.qual = qual
            self._input_sequence = self.sequence

        # sequence is already a Sequence object
        elif type(seq) == Sequence:
            self.id = id if id is not None else seq.id
            self.sequence = seq.sequence
            self.description = seq.description
            self.qual = seq.qual
            self._input_sequence = self.sequence
            self._annotations = seq._annotations

        # sequence is an iterable of (id, sequence, [qual])
        elif type(seq) in [list, tuple]:
            self.id = str(seq[0])
            self.sequence = str(seq[1]).upper()
            if len(seq) == 3 and qual is None:
                self.qual = seq[2]
            else:
                self.qual = qual
            self._input_sequence = self.sequence

        # sequence is a biopython SeqRecord
        elif type(seq) == SeqRecord:
            if qual is None:
                if "phred_quality" in seq.letter_annotations:
                    qual = "".join(
                        chr(q + 33) for q in seq.letter_annotations["phred_quality"]
                    )
                elif "solexa_quality" in seq.letter_annotations:
                    qual = "".join(
                        chr(q + 64) for q in seq.letter_annotations["solexa_quality"]
                    )
            self.id = id if id is not None else str(seq.id)
            self.description = str(seq.description)
            self.sequence = str(seq.seq).upper()
            self.qual = qual
            self._input_sequence = self.sequence

        # sequence is an annotation dict
        elif type(seq) in [dict, OrderedDict]:
            if self.id_key is None:
                id_options = ["seq_id", "sequence_id"]
                if any([k in seq for k in id_options]):
                    self.id_key = [k for k in id_options if k in seq][0]
            if self.seq_key is None:
                seq_options = ["vdj_nt", "sequence_nt", "sequence"]
                if any([k in seq for k in seq_options]):
                    self.seq_key = [k for k in seq_options if k in seq][0]
            self.id = id if id is not None else seq.get(self.id_key, None)
            self.sequence = seq.get(self.seq_key, None)
            self.qual = qual
            self._input_sequence = self.sequence
            for k in seq:
                if k not in self.annotations:
                    self.annotations[k] = seq[k]


def reverse_complement(
    sequence: Union[str, "Sequence"], in_place: bool = False
) -> Union[str, "Sequence"]:
    """
    Returns the reverse complement of a nucleotide sequence.

    Parameters
    ----------
    sequence : Union[str, Sequence]
        Nucleotide sequence to be reverse complemented.

    in_place : bool, default=False
        If ``True``, the input sequence will be reverse complemented in place.
        If ``False``, the reverse complemented sequence will be returned as a ``str``.


    Returns
    --------
    Union[str, Sequence]
        If ``in_place`` is ``False``, the reverse complemented sequence will be returned as a ``str``.
        If ``in_place`` is ``True``, the input sequence will be reverse complemented in place.

    """
    s = Sequence(sequence)
    rc_dict = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "Y": "R",
        "R": "Y",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
    rc = "".join([rc_dict.get(res, res) for res in s.sequence[::-1]])
    if in_place:
        s.sequence = rc
        return s
    return rc


def translate(
    sequence: Sequence,
    sequence_key: Optional[str] = None,
    frame: int = 1,
    allow_dots: bool = False,
) -> str:
    """
    Translates a nucleotide sequence.

    Parameters
    ----------
    sequence : Sequence
        ``Sequence`` object to be translated. Required.

    sequence_key : str, default=None
        Name of the annotation field containg the sequence to be translated.
        If not provided, ``sequence.sequence`` is used.

    frame : int, default=1
        Reading frame to translate. Default is ``1``.

    allow_dots : bool, default=False
        If ``True``, an all-dot codon ("...") will be translated as a single dot (".").
        Useful when translating IMGT-gapped sequences.Default is ``False``.


    Returns
    -------
    translated : str
        Translated sequence.

    """
    if sequence_key is not None:
        seq = Sequence(sequence[sequence_key])
    else:
        seq = Sequence(sequence)
    if seq is None:
        return None
    if frame not in range(1, 4):
        raise ValueError(f"Invalid frame: {frame}. Must be 1, 2 or 3.")
    start = frame - 1
    end = len(seq) - (len(seq[start:]) % 3)
    seq = seq[start:end]
    translated = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        if len(codon) != 3:
            continue
        if codon == "---":
            aa = "-"
        elif allow_dots and codon == "...":
            aa = "."
        else:
            aa = codon_lookup.get(codon, "X")
        translated += aa
    return translated


def codon_optimize(
    sequence,
    sequence_key: Optional[str] = None,
    id_key: Optional[str] = None,
    frame: Optional[int] = None,
    as_string: bool = False,
    in_place: bool = False,
) -> Union[str, "Sequence"]:
    """
    Optimizes the codons of a nucleotide or amino acid sequence.

    Parameters
    ----------
    sequence : Union[str, Sequence]
        Nucleotide or amino acid sequence to be optimized.

    sequence_key : str, default=None
        Name of the annotation field containg the sequence to be optimized.
        If not provided, ``sequence.sequence`` is used.

    id_key : str, default=None
        Name of the annotation field containg the sequence ID.
        If not provided, ``sequence.id`` is used.

    frame : int, default=1
        Reading frame to translate. Default is ``1``.

        .. note::
            ``frame`` is ignored if the input sequence is an amino acid sequence.

    as_string : bool, default=False
        If ``True``, the optimized sequence will be returned as a ``str``.
        If ``False``, the optimized sequence will be returned as a ``Sequence`` object.

    in_place : bool, default=False
        If ``True``, the input ``Sequence`` object will be returned with the optimized sequence populating
        the ``sequence.codon_optimized`` property. If ``False``, the optimized sequence will be returned
        as a ``str`` (if ``as_string`` is ``True``) or a ``Sequence`` object (if ``as_string`` is ``False``).

    Returns
    -------
    Union[str, Sequence]
        If ``in_place`` is ``True``, the input ``Sequence`` object will be returned with the optimized sequence populating
        the ``Sequence.sequence`` property.
        If ``as_string`` is ``True``, the optimized sequence will be returned as a ``str``.
        If ``as_string`` is ``False``, the optimized sequence will be returned as a new ``Sequence`` object.

    """
    # process input
    if not isinstance(sequence, Sequence):
        sequence = Sequence(sequence)
    if sequence_key is not None:
        seq = (
            sequence[sequence_key]
            if sequence[sequence_key] is not None
            else sequence.sequence
        )
    else:
        seq = sequence.sequence
    if id_key is not None:
        name = sequence[id_key] if id_key is not None else sequence.id
    else:
        name = sequence.id

    # get DNA sequence
    if not all([res.upper() in ["A", "C", "G", "T", "N", "-"] for res in sequence]):
        dna_seq = dc.reverse_translate(seq)
    else:
        frame = frame if frame is not None else 1
        if frame not in range(1, 4):
            raise ValueError(f"Invalid frame: {frame}. Must be 1, 2 or 3.")
        start = frame - 1
        dna_seq = seq[start:]

    # optimize codons
    problem = dc.DnaOptimizationProblem(
        sequence=dna_seq,
        constraints=[
            dc.EnforceTranslation(),
            dc.EnforceGCContent(maxi=0.56),
            dc.EnforceGCContent(maxi=0.64, window=100),
            dc.UniquifyAllKmers(10),
        ],
        objectives=[dc.CodonOptimize(species="h_sapiens")],
        logger=None,
    )
    problem.resolve_constraints(final_check=True)
    problem.optimize()
    if as_string:
        return str(problem.sequence)
    if in_place:
        sequence.sequence = str(problem.sequence)
        return sequence
    return Sequence(problem.sequence, id=name)


def read_json(
    json_file: str,
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "seq_id",
    sequence_key: str = "vdj_nt",
) -> Iterable[Sequence]:
    """
    Reads a JSON-formatted annotation file and returns ``Sequence`` objects.

    Parameters
    ----------
    json_file : str
        Path to the input JSON file. Required.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'locus'`` field is not ``'IGH'``:

        .. code-block:: python

            {'locus': 'IGH'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Sequence``
        objects. Fields must be column names in the input file.

    id_key : str, default="seq_id"
        Name of the annotation field containing the sequence ID. Used to
        populate the ``Sequence.id`` property.

    sequence_key : str, default="vdj_nt"
        Name of the annotation field containg the sequence. Used to
        populate the ``Sequence.sequence`` property.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    if match is None:
        match = {}
    if fields is not None:
        if id_key not in fields:
            fields.append(id_key)
        if sequence_key not in fields:
            fields.append(sequence_key)
    sequences = []
    with open(json_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            j = json.loads(line.strip())
            try:
                if all(
                    [nested_dict_lookup(j, k.split(".")) == v for k, v in match.items()]
                ):
                    if fields is not None:
                        _fields = [f for f in fields if f in j]
                        j = {f: j[f] for f in _fields}
                    sequences.append(Sequence(j, id_key=id_key, seq_key=sequence_key))
            except KeyError:
                continue
    return sequences


def read_csv(
    csv_file: str,
    separator: str = ",",
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Sequence]:
    """
    Reads a tabular annotation file and returns ``Sequence`` objects.

    Parameters
    ----------
    csv_file : str
        Path to the input CSV file. Required.

    separator : str, default=","
        Column separator. Default is ``","``.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'locus'`` field is not ``'IGH'``:

        .. code-block:: python

            {'locus': 'IGH'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Sequence``
        objects. Fields must be column names in the input file.

    id_key : str, default="sequence_id"
        Name of the annotation field containing the sequence ID. Used to
        populate the ``Sequence.id`` property.

    sequence_key : str, default="sequence"
        Name of the annotation field containg the sequence. Used to
        populate the ``Sequence.sequence`` property.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    if match is None:
        match = {}
    if fields is not None:
        if id_key not in fields:
            fields.append(id_key)
        if sequence_key not in fields:
            fields.append(sequence_key)
    sequences = []
    # df = pd.read_csv(csv_file, delimiter=delimiter)
    df = pl.read_csv(csv_file, separator=separator)
    # for _, r in df.iterrows():
    for r in df.iter_rows(named=True):
        # r = r.to_dict()
        try:
            if all([r[k] == v for k, v in match.items()]):
                if fields is not None:
                    _fields = [f for f in fields if f in r]
                    r = {f: r[f] for f in _fields}
                sequences.append(Sequence(r, id_key=id_key, seq_key=sequence_key))
        except KeyError:
            continue
    return sequences


def read_airr(
    tsv_file: str, match: Optional[dict] = None, fields: Optional[Iterable] = None
) -> Iterable[Sequence]:
    """
    Reads an AIRR-formatted annotation file and returns ``Sequence`` objects.

    Parameters
    ----------
    tsv_file : str
        Path to the input TSV file. Required.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'locus'`` field is not ``'IGH'``:

        .. code-block:: python

            {'locus': 'IGH'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Sequence``
        objects. Fields must be column names in the input file.


    Returns
    -------
    sequences : list of ``Sequences``
    """
    return read_csv(tsv_file, separator="\t", match=match, fields=fields)


def read_parquet(
    parquet_file: str,
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Sequence]:
    """
    Reads a Parquet file and returns ``Sequence`` objects.

    Parameters
    ----------
    parquet_file : str
        Path to the input Parquet file. Required.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'locus'`` field is not ``'IGH'``:

        .. code-block:: python

            {'locus': 'IGH'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Sequence``
        objects. Fields must be column names in the input file.

    id_key : str, default="sequence_id"
        Name of the annotation field containing the sequence ID. Used to
        populate the ``Sequence.id`` property.

    sequence_key : str, default="sequence"
        Name of the annotation field containg the sequence. Used to
        populate the ``Sequence.sequence`` property.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    if match is None:
        match = {}
    if fields is not None:
        if id_key not in fields:
            fields.append(id_key)
        if sequence_key not in fields:
            fields.append(sequence_key)
    sequences = []
    df = pl.read_parquet(parquet_file)
    for r in df.iter_rows(named=True):
        try:
            if all([r[k] == v for k, v in match.items()]):
                if fields is not None:
                    _fields = [f for f in fields if f in r]
                    r = {f: r[f] for f in _fields}
                sequences.append(Sequence(r, id_key=id_key, seq_key=sequence_key))
        except KeyError:
            continue
    return sequences


def sequences_from_polars(
    df: Union[pl.DataFrame, pl.LazyFrame],
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Sequence]:
    """
    Reads a Polars DataFrame and returns ``Sequence`` objects.

    Parameters
    ----------
    df : Union[pl.DataFrame, pl.LazyFrame]
        The input Polars DataFrame. Required.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'locus'`` field is not ``'IGH'``:

        .. code-block:: python

            {'locus': 'IGH'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Sequence``
        objects. Fields must be column names in the input file.

    id_key : str, default="sequence_id"
        Name of the annotation field containing the sequence ID. Used to
        populate the ``Sequence.id`` property.

    sequence_key : str, default="sequence"
        Name of the annotation field containg the sequence. Used to
        populate the ``Sequence.sequence`` property.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    if match is None:
        match = {}
    if fields is not None:
        if id_key not in fields:
            fields.append(id_key)
        if sequence_key not in fields:
            fields.append(sequence_key)
    sequences = []
    if isinstance(df, pl.LazyFrame):
        df = df.collect()
    for r in df.iter_rows(named=True):
        try:
            if all([r[k] == v for k, v in match.items()]):
                if fields is not None:
                    _fields = [f for f in fields if f in r]
                    r = {f: r[f] for f in _fields}
                sequences.append(Sequence(r, id_key=id_key, seq_key=sequence_key))
        except KeyError:
            continue
    return sequences


def sequences_from_pandas(
    df: pd.DataFrame,
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Sequence]:
    """
    Reads a pandas DataFrame and returns ``Sequence`` objects.

    Parameters
    ----------
    df : pd.DataFrame
        The input pandas DataFrame. Required.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'locus'`` field is not ``'IGH'``:

        .. code-block:: python

            {'locus': 'IGH'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Sequence``
        objects. Fields must be column names in the input file.

    id_key : str, default="sequence_id"
        Name of the annotation field containing the sequence ID. Used to
        populate the ``Sequence.id`` property.

    sequence_key : str, default="sequence"
        Name of the annotation field containg the sequence. Used to
        populate the ``Sequence.sequence`` property.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    if match is None:
        match = {}
    if fields is not None:
        if id_key not in fields:
            fields.append(id_key)
        if sequence_key not in fields:
            fields.append(sequence_key)
    sequences = []
    for _, r in df.iterrows():
        try:
            if all([r[k] == v for k, v in match.items()]):
                if fields is not None:
                    _fields = [f for f in fields if f in r]
                    r = {f: r[f] for f in _fields}
                sequences.append(Sequence(r, id_key=id_key, seq_key=sequence_key))
        except KeyError:
            continue
    return sequences


def determine_fastx_format(fastx_file: str) -> str:
    """
    Get the format of a FASTA or FASTQ file.

    Parameters
    ----------
    fastx_file : str
        The path to the FASTA or FASTQ file. Can be gzip-compressed.

    Returns
    -------
    str
        The file format -- either "fasta" or "fastq". If the initial non-whitespace
        character in the input file isn't ">" or "@", ``None`` is returned.
    """
    fmt = None
    # return None if the file doesn't exist
    if not os.path.exists(fastx_file):
        return fmt
    # get the appropriate open function
    if fastx_file.endswith(".gz"):
        open_fn = gzip.open
        mode = "rt"
    else:
        open_fn = open
        mode = "r"
    # check file contents to determine format
    try:
        with open_fn(fastx_file, mode) as f:
            line = next(f)
            while not line.strip():
                line = next(f)
            if line.lstrip().startswith(">"):
                fmt = "fasta"
            elif line.lstrip().startswith("@"):
                fmt = "fastq"
    except StopIteration:
        pass
    return fmt


def read_fastx(fastx: str) -> Iterable[Sequence]:
    """
    Reads FASTA or FASTQ-formatted sequence data and returns ``Sequence`` objects.
    Gzipped files are supported.

    Parameters
    ----------
    fastx : str
        Path to a FASTA or FASTQ-formatted file, optionally gzip-compressed. Required.

    Returns
    -------
    sequences : list of ``Sequences``
    """
    # if input is a FASTA-formatted string, write to a temp file
    # because pyfastx only accepts file paths, not file-like objects,
    # which means we can't use StringIO to do it all in-memory
    if os.path.isfile(fastx):
        tmp = None
    else:
        tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")
        tmp.write(fastx)
        fastx = tmp.name

    # get the appropriate reader class
    fmt = determine_fastx_format(fastx)
    if fmt is None:
        raise ValueError(f"Unsupported file format for {fastx}")
    reader = pyfastx.Fasta if fmt == "fasta" else pyfastx.Fastq

    # read sequences
    sequences = []
    for seq_tuple in reader(fastx, build_index=False):
        sequences.append(Sequence(seq_tuple))

    # cleanup
    if tmp is not None:
        from ..io import delete_files

        delete_files(tmp.name)
    return sequences


def parse_fastx(fastx: str) -> Sequence:
    """
    Parses FASTA or FASTQ-formatted sequence data and returns ``Sequence`` objects.
    This differs from ``read_fastx`` in that it yields sequences one at a time
    rather than reading all into memory and returning a list. This method is
    safe for extremely large files that are potentially too large to fit into memory.

    Parameters
    ----------
    fastx : str
        Path to a FASTA or FASTQ-formatted file, optionally gzip-compressed. Required.

    Yields
    -------
    sequences : ``Sequence``
    """
    # get the appropriate reader class
    fmt = determine_fastx_format(fastx)
    if fmt is None:
        raise ValueError(f"Unsupported file format for {fastx}")
    reader = pyfastx.Fasta if fmt == "fasta" else pyfastx.Fastq

    # parse sequences
    for seq_tuple in reader(fastx, build_index=False):
        yield Sequence(seq_tuple)


def read_fasta(fasta: str) -> Iterable[Sequence]:
    """
    Reads FASTA-formatted sequence data and returns ``Sequence`` objects.
    Gzipped files are supported.

    Parameters
    ----------
    fasta : str
        Either a FASTA-formatted string or the path to a FASTA file.
        FASTA files can be gzip-compressed. Required.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    # if input is a FASTA-formatted string, write to a temp file
    # because pyfastx only accepts file paths, not file-like objects,
    # which means we can't use StringIO to do it all in-memory
    if os.path.isfile(fasta):
        tmp = None
    else:
        tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")
        tmp.write(fasta)
        fasta = tmp.name
    sequences = []
    for name, sequence in pyfastx.Fasta(fasta, build_index=False, full_name=True):
        sequences.append(Sequence(sequence, id=name))

    # cleanup
    if tmp is not None:
        from ..io import delete_files

        delete_files(tmp.name)
    return sequences

    # if os.path.isfile(fasta):
    #     with open(fasta) as f:
    #         sequences = [Sequence(s) for s in SeqIO.parse(f, "fasta")]
    # else:
    #     from io import StringIO

    #     f = StringIO(fasta)
    #     sequences = [Sequence(s) for s in SeqIO.parse(f, "fasta")]
    # return sequences


def parse_fasta(fasta: str) -> Sequence:
    """
    Parses FASTA-formatted sequence data and returns ``Sequence`` objects. This
    differs from ``read_fasta`` in that it yields sequences one at a time
    rather than reading all into memory and returning a list. This method is
    safe for extremely large files that are potentially too large to fit into memory.

    Parameters
    ----------
    fasta : str
        Path to a FASTA-formatted file, optionally gzip-compressed. Required.

    Yields
    -------
    sequences : ``Sequence``
    """
    for name, sequence in pyfastx.Fasta(fasta, build_index=False, full_name=True):
        yield Sequence(sequence, id=name)


def read_fastq(fastq: str) -> Iterable[Sequence]:
    """
    Reads FASTQ-formatted sequence data and returns ``Sequence`` objects.
    Gzipped files are supported.

    Parameters
    ----------
    fastq : str
        Either a FASTQ-formatted string or the path to a FASTQ file.
        FASTQ files can be gzip-compressed. Required.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    # if input is a FASTA-formatted string, write to a temp file
    # because pyfastx only accepts file paths, not file-like objects,
    # which means we can't use StringIO to do it all in-memory
    if os.path.isfile(fastq):
        tmp = None
    else:
        tmp = tempfile.NamedTemporaryFile(delete=False, mode="w")
        tmp.write(fastq)
        fastq = tmp.name
    sequences = []
    for name, sequence, qual in pyfastx.Fastq(fastq, build_index=False):
        sequences.append(Sequence(sequence, id=name, qual=qual))

    # cleanup
    if tmp is not None:
        from ..io import delete_files

        delete_files(tmp.name)
    return sequences

    # if os.path.isfile(fastq):
    #     if fastq.endswith(".gz"):
    #         open_func = gzip.open
    #         mode = "rt"
    #         encoding = "utf-8"
    #     else:
    #         open_func = open
    #         mode = "r"
    #         encoding = None
    #     with open_func(fastq, mode=mode, encoding=encoding) as f:
    #         sequences = [Sequence(s) for s in SeqIO.parse(f, "fastq")]
    # else:
    #     from io import StringIO

    #     f = StringIO(fastq)
    #     sequences = [Sequence(s) for s in SeqIO.parse(f, "fastq")]
    # return sequences


def parse_fastq(fastq: str) -> Sequence:
    """
    Parses FASTQ-formatted sequence data and returns ``Sequence`` objects. This
    differs from ``read_fasta`` in that it yields sequences one at a time
    rather than reading all into memory and returning a list. This method is
    safe for extremely large files that are potentially too large to fit into memory.

    Parameters
    ----------
    fastq : str
        Path to a FASTQ-formatted file, optionally gzip-compressed. Required.

    Yields
    -------
    sequences : ``Sequence``
    """
    for name, sequence, qual in pyfastx.Fastq(fastq, build_index=False):
        yield Sequence(sequence, id=name, qual=qual)


def from_mongodb(db, collection, match=None, project=None, limit=None):
    if match is None:
        match = {}
    if project is None:
        project = {}
    if limit is not None:
        results = db[collection].find(match, project, limit=limit)
    else:
        results = db[collection].find(match, project)
    return [Sequence(r) for r in results]


# def to_fasta(
#     sequences: Union[str, Iterable[Sequence], Iterable[SeqRecord], Iterable[Iterable]],
#     fasta_file: Optional[str] = None,
#     # as_string: bool = False,
#     id_key: Optional[str] = None,
#     sequence_key: Optional[str] = None,
#     tempfile_dir: Optional[str] = None,
#     append_chain: bool = True,
# ) -> str:
#     """
#     Writes sequences to a FASTA-formatted file or returns a FASTA-formatted string.

#     Parameters
#     ----------
#     sequences : Iterable[Sequence]
#         Accepts any of the following:
#             1. list of abutils ``Sequence`` and/or ``Pair`` objects
#             2. FASTA/Q-formatted string
#             3. path to a FASTA/Q-formatted file
#             4. list of BioPython ``SeqRecord`` objects
#             5. list of lists/tuples, of the format ``[sequence_id, sequence]``
#         Required.

#         .. note::
#             Processing a list containing a mixture of ``Sequence`` and/or ``Pair`` objects is supported.

#     fasta_file : str, default=None
#         Path to the output FASTA file. If neither `fasta_file` nor `tempfile_dir`
#         are provided, a FASTA-formatted string will be returned.

#     id_key : str, default=None
#         Name of the annotation field containing the sequence ID. If not provided,
#         ``sequence.id`` is used.

#     sequence_key : str, default=None
#         Name of the annotation field containg the sequence. If not provided,
#         ``sequence.sequence`` is used.

#     tempfile_dir : str, optional
#         If `fasta_file` is not provided, directory into which the tempfile
#         should be created. If the directory does not exist, it will be
#         created.

#     append_chain : bool, default=True
#         If ``True``, the chain (heavy or light) will be appended to the sequence name:
#         ``>MySequence_heavy``.

#         .. note::
#             This option is ignored unless a list containing ``Pair`` objects is provided.


#     Returns
#     --------
#     fasta : str
#         Path to a FASTA file or a FASTA-formatted string

#     """
#     if isinstance(sequences, str):
#         # if sequences is already a FASTA/Q-formatted file
#         if os.path.isfile(sequences):
#             fasta_string = "\n".join([s.fasta for s in parse_fasta(sequences)])
#         # if sequences is a FASTA-formatted string
#         else:
#             fasta_string = sequences
#     # if sequences is a list of Sequences
#     elif all([isinstance(s, (Sequence, Pair)) for s in sequences]):
#         fasta_strings = []
#         for s in sequences:
#             if isinstance(s, Pair):
#                 fasta_strings.append(
#                     s.fasta(
#                         name_field=id_key,
#                         sequence_field=sequence_key,
#                         append_chain=append_chain,
#                     )
#                 )
#             else:
#                 fasta_strings.append(
#                     s.as_fasta(name_field=id_key, seq_field=sequence_key)
#                 )
#         fasta_string = "\n".join(fasta_strings)
#     # anything else..
#     else:
#         fasta_string = "\n".join(
#             [Sequence(s, id_key=id_key, seq_key=sequence_key).fasta for s in sequences]
#         )
#     # output
#     if fasta_file is not None:
#         make_dir(os.path.dirname(fasta_file))
#         with open(fasta_file, "w") as f:
#             f.write(fasta_string)
#         return fasta_file
#     elif tempfile_dir is not None:
#         make_dir(tempfile_dir)
#         ff = tempfile.NamedTemporaryFile(dir=tempfile_dir, delete=False)
#         ff.close()
#         return ff.name
#     return fasta_file


def to_fastq(
    sequences: Union[str, Iterable[Sequence], Iterable[SeqRecord], Iterable[Iterable]],
    fastq_file: Optional[str] = None,
    as_string: bool = False,
    id_key: Optional[str] = None,
    sequence_key: Optional[str] = None,
    tempfile_dir: str = "/tmp",
) -> str:
    """
    Writes sequences to a FASTQ-formatted file or returns a FASTQ-formatted string.

    Parameters
    ----------
    sequences : Iterable[Sequence]
        An iterable of any of the following:
            1. list of abutils ``Sequence`` objects
            2. FASTQ-formatted string
            3. path to a FASTQ-formatted file
            4. list of BioPython ``SeqRecord`` objects
            5. list of lists/tuples, of the format ``[sequence_id, sequence]``
        Required.

    fastq_file : str, default=None
        Path to the output FASTQ file. If not provided and `as_string` is ``False``,
        a file will be created using ``tempfile.NamedTemporaryFile()``.

    as_string : bool, default=False
        Return a FASTA-formatted string rather than writing to file.

    id_key : str, default=None
        Name of the annotation field containing the sequence ID. If not provided,
        ``sequence.id`` is used.

    sequence_key : str, default=None
        Name of the annotation field containg the sequence. If not provided,
        ``sequence.sequence`` is used.

    tempfile_dir : str, default="/tmp"
        If `fasta_file` is not provided, directory into which the tempfile
        should be created. If the directory does not exist, it will be
        created. Default is "/tmp".


    Returns
    --------
    fasta : str
        Path to a FASTA file or a FASTA-formatted string

    """
    if isinstance(sequences, str):
        # if sequences is already a FASTQ-formatted file
        if os.path.isfile(sequences):
            if as_string:
                with open(sequences, "r") as f:
                    fastq_string = f.read()
                return fastq_string
            return sequences
        # if sequences is a FASTQ-formatted string
        else:
            fastq_string = sequences
    # if sequences is a list of Sequences
    elif all([type(s) == Sequence for s in sequences]):
        if not all([s.qual is not None for s in sequences]):
            err = "Not all provided sequences have quality scores, which are required for FASTQ format."
            raise RuntimeError(err)
        ids = [s.get(id_key, s.id) if id_key is not None else s.id for s in sequences]
        seqs = [
            s.get(sequence_key, s.sequence) if sequence_key is not None else s.sequence
            for s in sequences
        ]
        quals = [s.qual for s in sequences]
        fastq_string = "\n".join(
            f"@{i}\n{s}\n+\n{q}" for i, s, q in zip(ids, seqs, quals)
        )
    # anything else..
    else:
        fastq_string = "\n".join(
            [Sequence(s, id_key=id_key, seq_key=sequence_key).fastq for s in sequences]
        )
    # output
    if as_string:
        return fastq_string
    if fastq_file is None:
        make_dir(tempfile_dir)
        ff = tempfile.NamedTemporaryFile(dir=tempfile_dir, delete=False)
        ff.close()
        fastq_file = ff.name
    else:
        make_dir(os.path.dirname(fastq_file))
    with open(fastq_file, "w") as f:
        f.write(fastq_string)
    return fastq_file


def sequences_to_polars(
    sequences: Iterable[Sequence],
    annotations: Optional[Iterable[str]] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> pl.DataFrame:
    """
    Converts a list of ``Sequence`` objects to a ``polars.DataFrame`` object.

    Parameters
    ----------
    sequences : Iterable[Sequence]
        List of ``Sequence`` objects to be converted to a ``polars.DataFrame`` object. Required.

    annotations : list, default=None
        A list of annotation fields to be included in the ``polars.DataFrame`` object.

    columns : list, default=None
        Deprecated, use ``annotations`` instead. If provided and ``annotations`` is not,
        ``annotations`` will be set to ``columns``. If both are provided, ``annotations``
        will supercede ``columns``.

    properties : list, default=None
        A list of properties to be included in the ``polars.DataFrame`` object.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the ``polars.DataFrame`` object.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the ``polars.DataFrame`` object.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the ``polars.DataFrame`` object.

    leading : str or list, default=None
        Field or list of fields to appear first in the ``polars.DataFrame`` object. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the ``polars.DataFrame`` object and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    pl.DataFrame
        A ``polars.DataFrame`` object.
    """
    data = []
    # read Sequence data
    for s in sequences:
        if not s.annotations:
            d = {"sequence_id": s.id, "sequence": s.sequence}
        else:
            d = s.annotations
        if properties is not None:
            for prop in properties:
                try:
                    d[prop] = getattr(s, prop)
                except AttributeError:
                    continue
        data.append(d)

    # populate DataFrame
    df = pl.DataFrame(data, infer_schema_length=None)

    # drop NaN
    if drop_na_columns:
        df = df[[s.name for s in df if not (s.null_count() == df.height)]]

    # excluded columns
    if exclude is not None:
        if isinstance(exclude, str):
            exclude = []
        cols = [c for c in df.columns if c not in exclude]
        df = df.select(cols)

    # reorder
    if order is not None:
        cols = [o for o in order if o in df.columns]
        df = df.select(cols)

    # leading columns
    if leading is not None:
        if isinstance(leading, str):
            leading = [leading]
        leading = [l for l in leading if l in df.columns]
        cols = leading + [c for c in df.columns if c not in leading]
        df = df.select(cols)

    if annotations is None and columns is not None:
        annotations = columns
    if annotations is not None:
        cols = [c for c in annotations if c in df.columns]
        df = df.select(cols)

    return df


def sequences_to_pandas(
    sequences: Iterable[Sequence],
    annotations: Optional[Iterable[str]] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> pd.DataFrame:
    """
    Converts a list of ``Sequence`` objects to a ``pandas.DataFrame`` object.

    Parameters
    ----------
    sequences : Iterable[Sequence]
        List of ``Sequence`` objects to be converted to a ``pandas.DataFrame`` object. Required.

    annotations : list, default=None
        A list of annotation fields to be included in the ``pandas.DataFrame`` object.

    columns : list, default=None
        Deprecated, use ``annotations`` instead. If provided and ``annotations`` is not,
        ``annotations`` will be set to ``columns``. If both are provided, ``annotations``
        will supercede ``columns``.

    properties : list, default=None
        A list of properties to be included in the ``pandas.DataFrame`` object.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the ``pandas.DataFrame`` object.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the `pandas.DataFrame`` object.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the ``pandas.DataFrame`` object.

    leading : str or list, default=None
        Field or list of fields to appear first in the ``pandas.DataFrame`` object. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the ``pandas.DataFrame`` object and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    pandas.DataFrame
        A ``pandas.DataFrame`` object.
    """
    df = sequences_to_polars(
        sequences,
        annotations=annotations,
        columns=columns,
        properties=properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )
    return df.to_pandas()


def sequences_to_csv(
    sequences: Iterable[Sequence],
    csv_file: str,
    separator: str = ",",
    header: bool = True,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Sequence`` objects to a CSV file.

    Parameters
    ----------
    sequences : Iterable[Sequence]
        List of ``Sequence`` objects to be saved to a CSV file. Required.

    csv_file : str
        Path to the output CSV file. Required.

    separator : str, default=","
        Column separator. Default is ``","``.

    header : bool, default=True
        If ``True``, the CSV file will contain a header row. Default is ``True``.

    columns : list, default=None
        A list of fields to be retained in the output CSV file. Fields must be column
        names in the input file.

    properties : list, default=None
        A list of properties to be included in the CSV file.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the CSV file.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the CSV file.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the CSV file.

    leading : str or list, default=None
        Field or list of fields to appear first in the CSV file. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the CSV file and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    None

    """
    df = sequences_to_polars(
        sequences,
        columns=columns,
        properties=properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )
    df.write_csv(csv_file, separator=separator, include_header=header)


def sequences_to_parquet(
    sequences: Iterable[Sequence],
    parquet_file: str,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Sequence`` objects to a Parquet file.

    Parameters
    ----------
    sequences : Iterable[Sequence]
        List of ``Sequence`` objects to be saved to a Parquet file. Required.

    parquet_file : str
        Path to the output Parquet file. Required.

    columns : list, default=None
        A list of fields to be retained in the output Parquet file. Fields must be column
        names in the input file.

    properties : list, default=None
    properties : list, default=None
        A list of properties to be included in the Parquet file.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the Parquet file.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the Parquet file.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the Parquet file.

    leading : str or list, default=None
        Field or list of fields to appear first in the Parquet file. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the Parquet file and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    None

    """
    df = sequences_to_polars(
        sequences,
        columns=columns,
        properties=properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )
    df.write_parquet(parquet_file)


def to_airr(
    sequences: Iterable[Sequence],
    airr_file: str,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Sequence`` objects to a CSV file.

    Parameters
    ----------
    sequences : Iterable[Sequence]
        List of ``Sequence`` objects to be saved to an AIRR file. Required.

    airr_file : str
        Path to the output AIRR file. Required.

    columns : list, default=None
        A list of fields to be retained in the output CSV file. Fields must be column
        names in the input file.

    properties : list, default=None
        A list of properties to be included in the CSV file.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the CSV file.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the CSV file.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the CSV file.

    leading : str or list, default=None
        Field or list of fields to appear first in the CSV file. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the CSV file and
        remaining fields will appear in the order provided in ``order``.

    Returns
    -------
    None

    """
    df = sequences_to_polars(
        sequences,
        columns=columns,
        properties=properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )
    df.write_csv(airr_file, separator="\t", include_header=True)
