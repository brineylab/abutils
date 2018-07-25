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
import json
import operator
import sys
import uuid

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..utils.utilities import nested_dict_lookup

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
    from functools import reduce
else:
    STR_TYPES = [str, unicode]


class Sequence(object):
    """
    Container for biological (RNA and DNA) sequences.

    ``seq`` can be one of several things:

        1) a raw sequence, as a string

        2) an iterable, formatted as ``[seq_id, sequence]``

        3) a dict, containing at least the ID (default key = 'seq_id') and a
           sequence (default key = 'vdj_nt'). Alternate ``id_key`` and ``seq_key``
           can be provided at instantiation.

        4) a Biopython ``SeqRecord`` object

        5) an AbTools ``Sequence`` object

    If ``seq`` is provided as a string, the sequence ID can optionally be
    provided via ``id``.  If ``seq`` is a string and ``id`` is not provided,
    a random sequence ID will be generated with ``uuid.uuid4()``.

    Quality scores can be supplied with ``qual`` or as part of a ``SeqRecord`` object.
    If providing both a SeqRecord object with quality scores and quality scores
    via ``qual``, the ``qual`` scores will override the SeqRecord quality scores.

    If ``seq`` is a dictionary, typically the result of a MongoDB query, the dictionary
    can be accessed directly from the ``Sequence`` instance. To retrive the value
    for ``'junc_aa'`` in the instantiating dictionary, you would simply::

        s = Sequence(dict)
        junc = s['junc_aa']

    If ``seq`` is a dictionary, an optional ``id_key`` and ``seq_key`` can be provided,
    which tells the ``Sequence`` object which field to use to populate ``Sequence.id`` and
    ``Sequence.sequence``. Defaults are ``id_key='seq_id'`` and ``seq_key='vdj_nt'``.

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
    def __init__(self, seq, id=None, qual=None, id_key='seq_id', seq_key='vdj_nt'):
        super(Sequence, self).__init__()
        self._input_sequence = None
        self._input_id = id
        self._input_qual = qual
        self.id = None
        self.sequence = None
        self.qual = None
        self._annotations = None
        self._fasta = None
        self._fastq = None
        self._reverse_complement = None
        self._strand = None
        self.id_key = id_key
        self.seq_key = seq_key

        self._process_input(seq, id, qual)


    def __len__(self):
        '''
        int: Returns the length of the ``sequence`` attribute
        '''
        return len(self.sequence)

    def __iter__(self):
        '''
        iter: Returns an iterator over the ``sequence`` attribute
        '''
        return iter(self.sequence)

    def __reversed__(self):
        '''
        str: Returns the reverse of the ``sequence`` attribute
        '''
        return ''.join(reversed(self.sequence))

    def __contains__(self, item):
        '''
        bool: If the instance was initialzed with a dictonary (which means
            the ``annotations`` attribute is not empty), ``__contains__(key)``
            will return ``True`` if ``key`` is in ``annotations.keys()``. If ``annotations``
            is an empty dict, indicating instantiation without a dictionary,
            ``__contains__(motif)`` will return True if ``motif`` is in the
            ``sequence attribute.

        '''
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
        if all([hasattr(other, 'sequence'), hasattr(other, 'id')]):
            return all([self.sequence == other.sequence, self.id == other.id])
        return False


    @property
    def fasta(self):
        '''
        str: Returns the sequence, as a FASTA-formatted string

        Note: The FASTA string is built using ``Sequence.id`` and ``Sequence.sequence``.
        '''
        if not self._fasta:
            self._fasta = '>{}\n{}'.format(self.id, self.sequence)
        return self._fasta

    @property
    def fastq(self):
        '''
        str: Returns the sequence, as a FASTQ-formatted string

        If ``Sequence.qual`` is ``None``, then ``None`` will be returned instead of a
        FASTQ string
        '''
        if self.qual is None:
            self._fastq = None
        else:
            if self._fastq is None:
                self._fastq = '@{}\n{}\n+\n{}'.format(self.id, self.sequence, self.qual)
        return self._fastq

    @property
    def reverse_complement(self):
        '''
        str: Returns the reverse complement of ``Sequence.sequence``.
        '''
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
            self._strand = 'plus'
        return self._strand

    @strand.setter
    def strand(self, strand):
        self._strand = strand


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
        return '>{}\n{}'.format(name, sequence)


    def region(self, start=0, end=None):
        '''
        Returns a region of ``Sequence.sequence``, in FASTA format.

        If called without kwargs, the entire sequence will be returned.

        Args:

            start (int): Start position of the region to be returned. Default
                is 0.

            end (int): End position of the region to be returned. Negative values
                will function as they do when slicing strings.

        Returns:

            str: A region of ``Sequence.sequence``, in FASTA format
        '''
        if end is None:
            end = len(self.sequence)
        return '>{}\n{}'.format(self.id, self.sequence[start:end])


    def keys(self):
        return self.annotations.keys()


    def values(self):
        return self.annotations.values()


    def get(self, key, default=None):
        return self.annotations.get(key, default)


    def _get_reverse_complement(self):
        rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
              'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W',
              'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
              'H': 'D', 'V': 'B', 'N': 'N'}
        return ''.join([rc.get(res, res) for res in self.sequence[::-1]])


    def _process_input(self, seq, id, qual):
        if type(seq) in STR_TYPES:
            if id is None:
                id = uuid.uuid4()
            if sys.version_info[0] == 2:
                self.id = id.encode('ascii') if isinstance(id, unicode) else id
                self.sequence = str(seq).encode('ascii').upper() if isinstance(str(seq), unicode) else str(seq).upper()
                self.qual = qual.encode('ascii') if isinstance(qual, unicode) else qual
            else:
                self.sequence = str(seq).upper()
                self.id = id
                self.qual = qual
            self._input_sequence = self.sequence
        elif type(seq) == Sequence:
            if sys.version_info[0] == 2:
                self.id = seq.id.encode('ascii') if isinstance(seq.id, unicode) else seq.id
                self.sequence = seq.sequence.encode('ascii') if isinstance(seq.sequence, unicode) else seq.sequence
                self.qual = seq.qual.encode('ascii') if isinstance(seq.qual, unicode) else seq.qual
            else:
                self.id = seq.id
                self.sequence = seq.sequence
                self.qual = seq.qual
            self._input_sequence = self.sequence
            self._annotations = seq._annotations
        elif type(seq) in [list, tuple]:
            if sys.version_info[0] == 2:
                self.id = str(seq[0]).encode('ascii').upper() if isinstance(str(seq[0]), unicode) else str(seq[0])
                self.sequence = str(seq[1]).encode('ascii').upper() if isinstance(str(seq[1]), unicode) else str(seq[1]).upper()
                self.qual = self.qual.encode('ascii') if isinstance(qual, unicode) else qual
            else:
                self.id = str(seq[0])
                self.sequence = str(seq[1]).upper()
                self.qual = qual
            self._input_sequence = self.sequence
        elif type(seq) == SeqRecord:
            if qual is None:
                if 'phred_quality' in seq.letter_annotations:
                    qual = seq.letter_annotations['phred_quality']
                elif 'solexa_quality' in seq.letter_annotations:
                    qual = seq.letter_annotations['solexa_quality']
            if sys.version_info[0] == 2:
                self.id = str(seq.id).encode('ascii') if isinstance(str(seq.id), unicode) else str(seq.id)
                self.sequence = str(seq.seq).encode('ascii').upper() if isinstance(str(seq.seq), unicode) else str(seq.seq).upper()
                self.qual = qual.encode('ascii') if isinstance(qual, unicode) else qual
            else:
                self.id = str(seq.id)
                self.sequence = str(seq.seq).upper()
                self.qual = qual
            self._input_sequence = self.sequence
        elif type(seq) in [dict, OrderedDict]:
            if sys.version_info[0] == 2:
                self.id = seq[self.id_key].encode('ascii') if isinstance(seq[self.id_key], unicode) else seq[self.id_key]
                self.sequence = seq[self.seq_key].encode('ascii').upper() if isinstance(seq[self.seq_key], unicode) else seq[self.seq_key].upper()
                self.qual = qual if isinstance(qual, unicode) else qual
            else:
                self.id = seq[self.id_key]
                self.sequence = seq[self.seq_key]
                self.qual = qual
            self._input_sequence = self.sequence
            self._annotations = seq


def read_json(json_file, match=None):
    if match is None:
        match = {}
    sequences = []
    with open(json_file, 'r') as f:
        for line in f:
            j = json.loads(line.strip())
            try:
                if all([nested_dict_lookup(j, k.split('.')) == v for k, v in match.items()]):
                    sequences.append(Sequence(j))
            except KeyError:
                continue
    return sequences


def read_fasta(fasta_file):
    with open(fasta_file) as f:
        sequences = [Sequence(s) for s in SeqIO.parse(f, 'fasta')]
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

            





