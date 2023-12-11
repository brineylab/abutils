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


from __future__ import absolute_import, division, print_function, unicode_literals

from collections import OrderedDict
import csv
import json
import operator
import os
import sys
import tempfile
from typing import Iterable, Optional, Union
import uuid

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pandas as pd

from ..utils.codons import codon_lookup
from ..utils.pipeline import make_dir
from ..utils.utilities import nested_dict_lookup

if sys.version_info[0] > 2:
    STR_TYPES = [
        str,
    ]
    from functools import reduce
else:
    STR_TYPES = [str, unicode]


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
    for ``'junc_aa'`` in the instantiating dictionary, you would simply::

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
    ``sequence`` attribute. An example of the distinction::

        d = {'name': 'MySequence', 'sequence': 'ATGC'}
        seq = Sequence(d, id_key='name', seq_key='sequence')

        seq['name']  # 'MySequence'
        seq[:2]  # 'AT'

    If the ``Sequence`` is instantiated with a dictionary, calls to ``__contains__()`` will
    return ``True`` if the supplied item is a key in the dictionary. In non-dict instantiations,
    ``__contains__()`` will look in the ``Sequence.sequence`` field directly (essentially a
    motif search). For example::

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
        seq: Union[str, Iterable, dict, SeqRecord],
        id: Optional[str] = None,
        qual: Optional[str] = None,
        annotations: Optional[dict] = None,
        id_key: Optional[str] = None,
        seq_key: Optional[str] = None,
    ):
        super(Sequence, self).__init__()
        self._input_sequence = None
        self._input_id = id
        self._input_qual = qual
        self._annotations = annotations
        self.id = None
        self.sequence = None
        self.qual = None
        self._fasta = None
        self._fastq = None
        self._reverse_complement = None
        self._strand = None
        self.id_key = id_key
        self.seq_key = seq_key
        self.description = None

        self._process_input(seq, id, qual)

    def __len__(self):
        """
        int: Returns the length of the ``sequence`` attribute
        """
        return len(self.sequence)

    def __iter__(self):
        """iter: Returns an iterator over the ``sequence`` attribute"""
        return iter(self.sequence)

    def __reversed__(self):
        """
        str: Returns the reverse of the ``sequence`` attribute
        """
        return "".join(reversed(self.sequence))

    def __contains__(self, item):
        """
        bool: If the instance was initialzed with a dictonary (which means
            the ``annotations`` attribute is not empty), ``__contains__(key)``
            will return ``True`` if ``key`` is in ``annotations.keys()``. If ``annotations``
            is an empty dict, indicating instantiation without a dictionary,
            ``__contains__(motif)`` will return True if ``motif`` is in the
            ``sequence attribute.

        """
        if self.annotations:
            return item in self.annotations.keys()
        return item in self.sequence

    def __getitem__(self, key):
        if type(key) == slice:
            return self.sequence[key]
        elif key in self.annotations.keys():
            return self.annotations.get(key, None)
        elif type(key) == int:
            return self.sequence[key]
        return None

    def __setitem__(self, key, val):
        self.annotations[key] = val

    def __eq__(self, other):
        if all([hasattr(other, "sequence"), hasattr(other, "id")]):
            return all([self.sequence == other.sequence, self.id == other.id])
        return False

    @classmethod
    def from_json(cls, data):
        if isinstance(data, str):
            data = json.loads(data)
        s = cls(None)
        s.__dict__ = data
        return s

    @property
    def fasta(self):
        """
        str: Returns the sequence, as a FASTA-formatted string

        Note: The FASTA string is built using ``Sequence.id`` and ``Sequence.sequence``.
        """
        if not self._fasta:
            self._fasta = ">{}\n{}".format(self.id, self.sequence)
        return self._fasta

    @property
    def fastq(self):
        """
        str: Returns the sequence, as a FASTQ-formatted string

        If ``Sequence.qual`` is ``None``, then ``None`` will be returned instead of a
        FASTQ string
        """
        if self.qual is None:
            self._fastq = None
        else:
            if self._fastq is None:
                self._fastq = "@{}\n{}\n+\n{}".format(self.id, self.sequence, self.qual)
        return self._fastq

    @property
    def reverse_complement(self):
        """
        str: Returns the reverse complement of ``Sequence.sequence``.
        """
        if self._reverse_complement is None:
            self._reverse_complement = self._get_reverse_complement()
        return self._reverse_complement

    @property
    def annotations(self):
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

    def translate(self, sequence_key=None, frame=1):
        return translate(self, sequence_key, frame)

    def as_fasta(self, name_field=None, seq_field=None):
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

        Args:

            start (int): Start position of the region to be returned. Default
                is 0.

            end (int): End position of the region to be returned. Negative values
                will function as they do when slicing strings.

        Returns:

            str: A region of ``Sequence.sequence``, in FASTA format
        """
        if end is None:
            end = len(self.sequence)
        return ">{}\n{}".format(self.id, self.sequence[start:end])

    def keys(self):
        return self.annotations.keys()

    def values(self):
        return self.annotations.values()

    def get(self, key, default=None):
        return self.annotations.get(key, default)

    def _get_reverse_complement(self):
        rc = {
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
        return "".join([rc.get(res, res) for res in self.sequence[::-1]])

    def _process_input(self, seq, id, qual):
        if type(seq) in STR_TYPES:
            if id is None:
                id = uuid.uuid4()
            self.sequence = str(seq).upper()
            self.id = id
            self.qual = qual
            self._input_sequence = self.sequence
        elif type(seq) == Sequence:
            self.id = id if id is not None else seq.id
            self.sequence = seq.sequence
            self.description = seq.description
            self.qual = seq.qual
            self._input_sequence = self.sequence
            self._annotations = seq._annotations
        elif type(seq) in [list, tuple]:
            self.id = str(seq[0])
            self.sequence = str(seq[1]).upper()
            self.qual = qual
            self._input_sequence = self.sequence
        elif type(seq) == SeqRecord:
            if qual is None:
                if "phred_quality" in seq.letter_annotations:
                    qual = "".join(chr(q + 33) for q in seq.letter_annotations["phred_quality"])
                elif "solexa_quality" in seq.letter_annotations:
                    qual = "".join(chr(q + 64) for q in seq.letter_annotations["solexa_quality"])
            self.id = id if id is not None else str(seq.id)
            self.description = str(seq.description)
            self.sequence = str(seq.seq).upper()
            self.qual = qual
            self._input_sequence = self.sequence
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


def translate(
    sequence: Sequence, sequence_key: Optional[str] = None, frame: int = 1
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


    Returns
    -------
    translated : str
    """
    if sequence_key is not None:
        seq = Sequence(
            nested_dict_lookup(sequence.annotations, sequence_key.split("."))
        )
    else:
        seq = Sequence(sequence)
    if seq is None:
        return None
    start = (frame % 3) - 1
    end = len(seq) - (len(seq[start:]) % 3)
    seq = seq[start:end]
    translated = ""
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        if all([c == "-" for c in codon]):
            aa = "-"
        else:
            aa = codon_lookup.get(codon, "X")
        translated += aa
    return translated


def read_json(
    json_file: str,
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "seq_id",
    sequence_key: str = "vdj_nt",
):
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
        the ``'locus'`` field is no ``'IGH'``:
            ``{'locus': 'IGH'}

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
    delimiter: str = ",",
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

    delimiter : str, default=","
        Column delimiter. Default is ``","``.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'locus'`` field is no ``'IGH'``:
            ``{'locus': 'IGH'}

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
    df = pd.read_csv(csv_file, delimiter=delimiter)
    for _, r in df.iterrows():
        r = r.to_dict()
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
        the ``'locus'`` field is no ``'IGH'``:
            ``{'locus': 'IGH'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Sequence``
        objects. Fields must be column names in the input file.


    Returns
    -------
    sequences : list of ``Sequences``
    """
    return read_csv(tsv_file, delimiter="\t", match=match, fields=fields)


def read_fasta(fasta_file: str) -> Iterable[Sequence]:
    """
    Reads a FASTA-formatted  file and returns ``Sequence`` objects.

    Parameters
    ----------
    fasta_file : str
        Path to the input FASTA file. Required.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    with open(fasta_file) as f:
        sequences = [Sequence(s) for s in SeqIO.parse(f, "fasta")]
    return sequences


def read_fastq(fastq_file):
    """
    Reads a FASTQ-formatted  file and returns ``Sequence`` objects. 
    Gzipped files are supported.

    Parameters
    ----------
    fastq_file : str
        Path to the input FASTQ file. Required.


    Returns
    -------
    sequences : list of ``Sequences``

    """
    if fastq_file.endswith(".gz"):
        import gzip
        open_func = gzip.open
    else:
        open_func = open
    with open_func(fastq_file, "rt") as f:
        sequences = [Sequence(s) for s in SeqIO.parse(f, "fastq")]
    return sequences


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


def to_fasta(
    sequences: Union[str, Iterable[Sequence], Iterable[SeqRecord], Iterable[Iterable]],
    fasta_file: Optional[str] = None,
    as_string: bool = False,
    id_key: Optional[str] = None,
    sequence_key: Optional[str] = None,
    tempfile_dir: str = "/tmp",
) -> str:
    """
    Writes sequences to a FASTA-formatted file or returns a FASTA-formatted string.

    Parameters
    ----------
    sequences : Iterable[Sequence]
        An iterable of any of the following:
            1. list of abutils ``Sequence`` objects
            2. FASTA-formatted string
            3. path to a FASTA-formatted file
            4. list of BioPython ``SeqRecord`` objects
            5. list of lists/tuples, of the format ``[sequence_id, sequence]``
        Required.

    fasta_file : str, default=None
        Path to the output FASTA file. If not provided and `as_string` is ``False``,
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
        # if sequences is already a FASTA-formatted file
        if os.path.isfile(sequences):
            if as_string:
                with open(sequences, "r") as f:
                    fasta_string = f.read()
                return fasta_string
            return sequences
        # if sequences is a FASTA-formatted string
        else:
            fasta_string = sequences
    # if sequences is a list of Sequences
    elif all([type(s) == Sequence for s in sequences]):
        ids = [s.get(id_key, s.id) if id_key is not None else s.id for s in sequences]
        seqs = [
            s.get(sequence_key, s.sequence) if sequence_key is not None else s.sequence
            for s in sequences
        ]
        fasta_string = "\n".join(f">{i}\n{s}" for i, s in zip(ids, seqs))
    # anything else..
    else:
        fasta_string = "\n".join(
            [Sequence(s, id_key=id_key, seq_key=sequence_key).fasta for s in sequences]
        )
    # output
    if as_string:
        return fasta_string
    if fasta_file is None:
        make_dir(tempfile_dir)
        ff = tempfile.NamedTemporaryFile(dir=tempfile_dir, delete=False)
        ff.close()
        fasta_file = ff.name
    else:
        make_dir(os.path.dirname(fasta_file))
    with open(fasta_file, "w") as f:
        f.write(fasta_string)
    return fasta_file


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
        fastq_string = "\n".join(f"@{i}\n{s}\n+\n{q}" for i, s, q in zip(ids, seqs, quals))
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