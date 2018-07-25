#!/usr/bin/env python
# filename: pair.py


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

import copy
import sys
import traceback

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from .sequence import Sequence
from ..utils import germlines
from ..utils.alignment import global_alignment


if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
else:
    STR_TYPES = [str, unicode]


class Pair(object):
    '''
    Holds a pair of sequences, corresponding to HC and LC of a single mAb.

    Input is a list of dicts, with each dict containing sequence information from a single
    chain, formatted as would be returned from a query on a MongoDB database containing
    AbStar output.
    '''
    def __init__(self, seqs, name=None, h_selection_func=None, l_selection_func=None):
        self._seqs = seqs
        self._heavy = None
        self._light = None
        self._heavies = [s for s in seqs if s['chain'] == 'heavy']
        self._lights = [s for s in seqs if s['chain'] in ['kappa', 'lambda']]
        self._name = name
        self._fasta = None
        self._sample = None
        self._subject = None
        self._group = None
        self._experiment = None
        self._timepoint = None
        self._is_pair = None
        self._vrc01_like = None
        self._lineage = None
        self._select_heavy = h_selection_func
        self._select_light = l_selection_func

    def __eq__(self, other):
        return (self.heavy, self.light) == (other.heavy, other.light)

    def __ne__(self, other):
        return not self == other

    def __hash(self):
        return hash((self.heavy, self.light))


    @property
    def heavy(self):
        if self._heavy is None:
            if len(self._heavies) > 0:
                if self._select_heavy is not None:
                    h = self._select_heavy(self._heavies)
                    if all([h is not None, type(h) != Sequence]):
                        h = Sequence(h)
                    self._heavy = h
                else:
                    h = self._heavies[0]
                    if type(h) != Sequence:
                        h = Sequence(h)
                    self._heavy = h
            else:
                self._heavy = None
        return self._heavy

    @heavy.setter
    def heavy(self, heavy):
        self._heavy = heavy

    @property
    def light(self):
        if self._light is None:
            if len(self._lights) > 0:
                if self._select_light is not None:
                    l = self._select_light(self._lights)
                    if all([l is not None, type(l) != Sequence]):
                        l = Sequence(l)
                    self._light = l
                else:
                    l = self._lights[0]
                    if type(l) != Sequence:
                        l = Sequence(l)
                    self._light = l
            else:
                self._light = None
        return self._light

    @light.setter
    def light(self, light):
        self._light = light

    @property
    def is_pair(self):
        if all([self.heavy is not None, self.light is not None]):
            return True
        return False

    @property
    def lineage(self):
        if self._lineage is None:
            self._lineage = self.heavy['clonify']['id']
        return self._lineage

    @property
    def vrc01_like(self):
        if self._vrc01_like is None:
            if any([self.heavy is None, self.light is None]):
                self._vrc01_like = False
            else:
                self._vrc01_like = all([self.heavy['v_gene']['gene'] == 'IGHV1-2', self.light['cdr3_len'] == 5])
        return self._vrc01_like

    @property
    def name(self):
        if self._name is None:
            if self.heavy is not None:
                self._name = self.heavy['seq_id']
            elif self.light is not None:
                self._name = self.light['seq_id']
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def sample(self):
        if self._sample is None:
            slist = []
            if self.experiment is not None:
                slist.append(str(self.experiment))
            if self.group is not None:
                slist.append(str(self.group))
            if self.subject is not None:
                slist.append(str(self.subject))
            if self.timepoint is not None:
                slist.append(str(self.timepoint))
            if slist:
                self._sample = '|'.join(slist)
        return self._sample

    @property
    def subject(self):
        if self._subject is None:
            if self.heavy is not None and 'subject' in self.heavy.keys():
                self._subject = self.heavy['subject']
            elif self.light is not None and 'subject' in self.light.keys():
                self._subject = self.light['subject']
        return self._subject

    @subject.setter
    def subject(self, subject):
        self._subject = subject

    @property
    def group(self):
        if self._group is None:
            if self.heavy is not None and 'group' in self.heavy.keys():
                self._group = self.heavy['group']
            elif self.light is not None and 'group' in self.light.keys():
                self._group = self.light['group']
        return self._group

    @group.setter
    def group(self, group):
        self._group = group

    @property
    def experiment(self):
        if self._experiment is None:
            if self.heavy is not None and 'experiment' in self.heavy.keys():
                self._experiment = self.heavy['experiment']
            elif self.light is not None and 'experiment' in self.light.keys():
                self._experiment = self.light['experiment']
        return self._experiment

    @experiment.setter
    def experiment(self, experiment):
        self._experiment = experiment

    @property
    def timepoint(self):
        if self._timepoint is None:
            if self.heavy is not None and 'timepoint' in self.heavy.keys():
                self._timepoint = self.heavy['timepoint']
            elif self.light is not None and 'timepoint' in self.light.keys():
                self._timepoint = self.light['timepoint']
        return self._timepoint

    @timepoint.setter
    def timepoint(self, timepoint):
        self._timepoint = timepoint


    def refine(self, heavy=True, light=True, species='human'):
        for seq in [s for s in [self.heavy, self.light] if s is not None]:
            try:
                self.remove_ambigs(seq)
                self._refine_v(seq, species)
                self._refine_j(seq, species)
                self._retranslate(seq)
            except:
                print('REFINEMENT FAILED: {}, {} chain'.format(seq['seq_id'], seq['chain']))
                print(traceback.format_exception_only(*sys.exc_info()[:2]))


    @staticmethod
    def remove_ambigs(seq):
        # fix Ns in the nucleotide sequence
        vdj = ''
        for s, g in zip(seq['vdj_nt'], seq['vdj_germ_nt']):
            if s.upper() == 'N':
                vdj += g
            else:
                vdj += s
        seq['vdj_nt'] = vdj
        # fix Xs in the amino acid sequence
        vdj = ''
        for s, g in zip(seq['vdj_aa'], seq['vdj_germ_aa']):
            if s.upper() == 'X':
                vdj += g
            else:
                vdj += s
        seq['vdj_aa'] = vdj

    @staticmethod
    def _refine_v(seq, species):
        '''
        Completes the 5' end of a a truncated sequence with germline nucleotides.
        Input is a MongoDB dict (seq) and the species.
        '''
        vgerm = germlines.get_germline(seq['v_gene']['full'], species)
        aln = global_alignment(seq['vdj_nt'], vgerm)
        prepend = ''
        for s, g in zip(aln.aligned_query, aln.aligned_target):
            if s != '-':
                break
            else:
                prepend += g
        seq['vdj_nt'] = prepend + seq['vdj_nt']

    @staticmethod
    def _refine_j(seq, species):
        '''
        Completes the 3' end of a a truncated sequence with germline nucleotides.
        Input is a MongoDB dict (seq) and the species.
        '''
        jgerm = germlines.get_germline(seq['j_gene']['full'], species)
        aln = global_alignment(seq['vdj_nt'], jgerm)
        append = ''
        for s, g in zip(aln.aligned_query[::-1], aln.aligned_target[::-1]):
            if s != '-':
                break
            else:
                append += g
        seq['vdj_nt'] = seq['vdj_nt'] + append[::-1]

    @staticmethod
    def _retranslate(seq):
        '''
        Retranslates a nucleotide sequence following refinement.
        Input is a Pair sequence (basically a dict of MongoDB output).
        '''
        if len(seq['vdj_nt']) % 3 != 0:
            trunc = len(seq['vdj_nt']) % 3
            seq['vdj_nt'] = seq['vdj_nt'][:-trunc]
        seq['vdj_aa'] = Seq(seq['vdj_nt'], generic_dna).translate()


    def fasta(self, key='vdj_nt', append_chain=True):
        '''
        Returns the sequence pair as a fasta string. If the Pair object contains
        both heavy and light chain sequences, both will be returned as a single string.

        By default, the fasta string contains the 'vdj_nt' sequence for each chain. To change,
        use the <key> option to select an alternate sequence.

        By default, the chain (heavy or light) will be appended to the sequence name:

        >MySequence_heavy

        To just use the pair name (which will result in duplicate sequence names for Pair objects
        with both heavy and light chains), set <append_chain> to False.
        '''
        fastas = []
        for s, chain in [(self.heavy, 'heavy'), (self.light, 'light')]:
            if s is not None:
                c = '_{}'.format(chain) if append_chain else ''
                fastas.append('>{}{}\n{}'.format(s['seq_id'], c, s[key]))
        return '\n'.join(fastas)


def get_pairs(db, collection, experiment=None, subject=None, group=None, name='seq_id',
    delim=None, delim_occurance=1, pairs_only=False, h_selection_func=None, l_selection_func=None):
    '''
    Gets sequences and assigns them to the appropriate mAb pair, based on the sequence name.

    Inputs:

    ::db:: is a pymongo database connection object
    ::collection:: is the collection name, as a string
    If ::subject:: is provided, only sequences with a 'subject' field matching ::subject:: will
        be included. ::subject:: can be either a single subject (as a string) or an iterable
        (list or tuple) of subject strings.
    If ::group:: is provided, only sequences with a 'group' field matching ::group:: will
        be included. ::group:: can be either a single group (as a string) or an iterable
        (list or tuple) of group strings.
    ::name:: is the dict key of the field to be used to group the sequences into pairs.
        Default is 'seq_id'
    ::delim:: is an optional delimiter used to truncate the contents of the ::name:: field.
        Default is None, which results in no name truncation.
    ::delim_occurance:: is the occurance of the delimiter at which to trim. Trimming is performed
        as delim.join(name.split(delim)[:delim_occurance]), so setting delim_occurance to -1 will
        trucate after the last occurance of delim. Default is 1.
    ::pairs_only:: setting to True results in only truly paired sequences (pair.is_pair == True)
        will be returned. Default is False.

    Returns a list of Pair objects, one for each mAb pair.
    '''
    match = {}
    if subject is not None:
        if type(subject) in (list, tuple):
            match['subject'] = {'$in': subject}
        elif type(subject) in STR_TYPES:
            match['subject'] = subject
    if group is not None:
        if type(group) in (list, tuple):
            match['group'] = {'$in': group}
        elif type(group) in STR_TYPES:
            match['group'] = group
    if experiment is not None:
        if type(experiment) in (list, tuple):
            match['experiment'] = {'$in': experiment}
        elif type(experiment) in STR_TYPES:
            match['experiment'] = experiment
    seqs = list(db[collection].find(match))
    return assign_pairs(seqs, name=name, delim=delim,
        delim_occurance=delim_occurance, pairs_only=pairs_only,
        h_selection_func=h_selection_func, l_selection_func=l_selection_func)


def assign_pairs(seqs, name='seq_id', delim=None, delim_occurance=1, pairs_only=False,
                 h_selection_func=None, l_selection_func=None):
    '''
    Assigns sequences to the appropriate mAb pair, based on the sequence name.

    Inputs:

    ::seqs:: is a list of dicts, of the format returned by querying a MongoDB containing
        Abstar output.
    ::name:: is the dict key of the field to be used to group the sequences into pairs.
        Default is 'seq_id'
    ::delim:: is an optional delimiter used to truncate the contents of the ::name:: field.
        Default is None, which results in no name truncation.
    ::delim_occurance:: is the occurance of the delimiter at which to trim. Trimming is performed
        as delim.join(name.split(delim)[:delim_occurance]), so setting delim_occurance to -1 will
        trucate after the last occurance of delim. Default is 1.
    ::pairs_only:: setting to True results in only truly paired sequences (pair.is_pair == True)
        will be returned. Default is False.

    Returns a list of Pair objects, one for each mAb pair.
    '''
    pdict = {}
    for s in seqs:
        if delim is not None:
            pname = delim.join(s[name].split(delim)[:delim_occurance])
        else:
            pname = s[name]
        if pname not in pdict:
            pdict[pname] = [s, ]
        else:
            pdict[pname].append(s)
    pairs = [Pair(pdict[n], name=n,
                  h_selection_func=h_selection_func,
                  l_selection_func=l_selection_func) for n in pdict.keys()]
    if pairs_only:
        pairs = [p for p in pairs if p.is_pair]
    return pairs


def deduplicate(pairs, aa=False, ignore_primer_regions=False):
    '''
    Removes duplicate sequences from a list of Pair objects.

    If a Pair has heavy and light chains, both chains must identically match heavy and light chains
    from another Pair to be considered a duplicate. If a Pair has only a single chain,
    identical matches to that chain will cause the single chain Pair to be considered a duplicate,
    even if the comparison Pair has both chains.

    Note that identical sequences are identified by simple string comparison, so sequences of
    different length that are identical over the entirety of the shorter sequence are not
    considered duplicates.

    By default, comparison is made on the nucleotide sequence. To use the amino acid sequence instead,
    set aa=True.
    '''
    nr_pairs = []
    just_pairs = [p for p in pairs if p.is_pair]
    single_chains = [p for p in pairs if not p.is_pair]
    _pairs = just_pairs + single_chains
    for p in _pairs:
        duplicates = []
        for nr in nr_pairs:
            identical = True
            vdj = 'vdj_aa' if aa else 'vdj_nt'
            offset = 4 if aa else 12
            if p.heavy is not None:
                if nr.heavy is None:
                    identical = False
                else:
                    heavy = p.heavy[vdj][offset:-offset] if ignore_primer_regions else p.heavy[vdj]
                    nr_heavy = nr.heavy[vdj][offset:-offset] if ignore_primer_regions else nr.heavy[vdj]
                    if heavy != nr_heavy:
                        identical = False
            if p.light is not None:
                if nr.light is None:
                    identical = False
                else:
                    light = p.light[vdj][offset:-offset] if ignore_primer_regions else p.light[vdj]
                    nr_light = nr.light[vdj][offset:-offset] if ignore_primer_regions else nr.light[vdj]
                    if light != nr_light:
                        identical = False
            duplicates.append(identical)
        if any(duplicates):
            continue
        else:
            nr_pairs.append(p)
    return nr_pairs


def refine(pairs, heavy=True, light=True, species='human'):
    refined_pairs = copy.deepcopy(pairs)
    for p in refined_pairs:
        p.refine(heavy, light, species)
    return refined_pairs
