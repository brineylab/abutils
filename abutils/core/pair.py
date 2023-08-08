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
import csv
import sys
import traceback

from Bio.Seq import Seq

from .sequence import Sequence
from ..utils import germlines
from ..utils.alignment import global_alignment
from ..utils.utilities import nested_dict_lookup


if sys.version_info[0] > 2:
    STR_TYPES = [
        str,
    ]
else:
    STR_TYPES = [str, unicode]


class Pair(object):
    """
    Holds a pair of sequences, corresponding to HC and LC of a single mAb.

    Input is a list of dicts, with each dict containing sequence information from a single
    chain, formatted as would be returned from a query on a MongoDB database containing
    AbStar output.
    """

    def __init__(self, seqs, name=None, chain_selection_func=None):
        self._seqs = seqs
        self._receptor = None
        self._heavy = None
        self._light = None
        self._alpha = None
        self._beta = None
        self._delta = None
        self._gamma = None
        self._heavies = None
        self._lights = None
        self._alphas = None
        self._betas = None
        self._deltas = None
        self._gammas = None
        self._name = name
        self._fasta = None
        self._sample = None
        self._subject = None
        self._group = None
        self._experiment = None
        self._timepoint = None
        self._is_pair = None
        self._lineage = None
        self._select_chain = (
            chain_selection_func
            if chain_selection_func is not None
            else self._chain_selector
        )

    def __eq__(self, other):
        return (self.heavy, self.light) == (other.heavy, other.light)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        name = self.name
        heavy = self.heavy.sequence if self.heavy is not None else None
        light = self.light.sequence if self.light is not None else None
        return hash((name, heavy, light))

    @property
    def receptor(self):
        if self._receptor is None:
            if all([s["chain"] is not None for s in self._seqs]):
                if all(
                    [s["chain"] in ["heavy", "kappa", "lambda"] for s in self._seqs]
                ):
                    self._receptor = "bcr"
                elif all(
                    [
                        s["chain"] in ["alpha", "beta", "delta", "gamma"]
                        for s in self._seqs
                    ]
                ):
                    self._receptor = "tcr"
                else:
                    self._receptor = "unknown"
            elif all([s["locus"] is not None for s in self._seqs]):
                if all(
                    [s["locus"].lower() in ["igh", "igk", "igl"] for s in self._seqs]
                ):
                    self._receptor = "bcr"
                elif all(
                    [
                        s["locus"].lower() in ["tra", "trb", "trd", "trg"]
                        for s in self._seqs
                    ]
                ):
                    self._receptor = "tcr"
                else:
                    self._receptor = "unknown"
        return self._receptor

    @receptor.setter
    def receptor(self, receptor):
        self._receptor = receptor

    @property
    def heavy(self):
        if self._heavy is None:
            self._heavy = self._select_chain(self.heavies)
        return self._heavy

    @heavy.setter
    def heavy(self, heavy):
        self._heavy = heavy

    @property
    def light(self):
        if self._light is None:
            self._light = self._select_chain(self.lights)
        return self._light

    @light.setter
    def light(self, light):
        self._light = light

    @property
    def alpha(self):
        if self._alpha is None:
            self._alpha = self._select_chain(self.alphas)
        return self._alpha

    @alpha.setter
    def alpha(self, alpha):
        self._alpha = alpha

    @property
    def beta(self):
        if self._beta is None:
            self._beta = self._select_chain(self.betas)
        return self._beta

    @beta.setter
    def beta(self, beta):
        self._beta = beta

    @property
    def delta(self):
        if self._delta is None:
            self._delta = self._select_chain(self.deltas)
        return self._delta

    @delta.setter
    def delta(self, delta):
        self._delta = delta

    @property
    def gamma(self):
        if self._gamma is None:
            self._gamma = self._select_chain(self.gammas)
        return self._gamma

    @gamma.setter
    def gamma(self, gamma):
        self._gamma = gamma

    @property
    def heavies(self):
        if self._heavies is None:
            if all([s["chain"] is not None for s in self._seqs]):
                self._heavies = [s for s in self._seqs if s["chain"] == "heavy"]
            elif all([s["locus"] is not None for s in self._seqs]):
                self._heavies = [s for s in self._seqs if s["locus"] == "IGH"]
        return self._heavies

    @property
    def lights(self):
        if self._lights is None:
            if all([s["chain"] is not None for s in self._seqs]):
                self._lights = [
                    s for s in self._seqs if s["chain"] in ["kappa", "lambda"]
                ]
            elif all([s["locus"] is not None for s in self._seqs]):
                self._lights = [s for s in self._seqs if s["locus"] in ["IGK", "IGL"]]
        return self._lights

    @property
    def alphas(self):
        if self._alphas is None:
            if all([s["chain"] is not None for s in self._seqs]):
                self._alphas = [s for s in self._seqs if s["chain"] == "alpha"]
            elif all([s["locus"] is not None for s in self._seqs]):
                self._alphas = [s for s in self._seqs if s["locus"] == "TRA"]
        return self._alphas

    @property
    def betas(self):
        if self._betas is None:
            if all([s["chain"] is not None for s in self._seqs]):
                self._betas = [s for s in self._seqs if s["chain"] == "beta"]
            elif all([s["locus"] is not None for s in self._seqs]):
                self._betas = [s for s in self._seqs if s["locus"] == "TRB"]
        return self._betas

    @property
    def deltas(self):
        if self._deltas is None:
            if all([s["chain"] is not None for s in self._seqs]):
                self._deltas = [s for s in self._seqs if s["chain"] == "delta"]
            elif all([s["locus"] is not None for s in self._seqs]):
                self._deltas = [s for s in self._seqs if s["locus"] == "TRD"]
        return self._deltas

    @property
    def gammas(self):
        if self._gammas is None:
            if all([s["chain"] is not None for s in self._seqs]):
                self._gammas = [s for s in self._seqs if s["chain"] == "gamma"]
            elif all([s["locus"] is not None for s in self._seqs]):
                self._gammas = [s for s in self._seqs if s["locus"] == "TRG"]
        return self._gammas

    @property
    def is_pair(self):
        if all([self.heavy is not None, self.light is not None]):
            return True
        elif all([self.alpha is not None, self.beta is not None]):
            return True
        elif all([self.gamma is not None, self.delta is not None]):
            return True
        return False

    @property
    def lineage(self):
        if self._lineage is None:
            if "clonify" in self.heavy:
                self._lineage = self.heavy["clonify"]["id"]
        return self._lineage

    @property
    def name(self):
        if self._name is None:
            if self.heavy is not None:
                if self.heavy.id is not None:
                    self._name = self.heavy.id
                elif self.heavy["sequence_id"] is not None:
                    self._name = self.heavy["sequence_id"]
                elif self.heavy["seq_id"] is not None:
                    self._name = self.heavy["seq_id"]
                else:
                    pass
            elif self.light is not None:
                if self.light.id is not None:
                    self._name = self.light.id
                elif self.light["sequence_id"] is not None:
                    self._name = self.light["sequence_id"]
                elif self.light["seq_id"] is not None:
                    self._name = self.light["seq_id"]
                else:
                    pass
            elif self.alpha is not None:
                if self.alpha.id is not None:
                    self._name = self.alpha.id
                elif self.alpha["sequence_id"] is not None:
                    self._name = self.alpha["sequence_id"]
                elif self.alpha["seq_id"] is not None:
                    self._name = self.alpha["seq_id"]
                else:
                    pass
            elif self.beta is not None:
                if self.beta.id is not None:
                    self._name = self.beta.id
                elif self.beta["sequence_id"] is not None:
                    self._name = self.beta["sequence_id"]
                elif self.beta["seq_id"] is not None:
                    self._name = self.beta["seq_id"]
                else:
                    pass
            elif self.delta is not None:
                if self.delta.id is not None:
                    self._name = self.delta.id
                elif self.delta["sequence_id"] is not None:
                    self._name = self.delta["sequence_id"]
                elif self.delta["seq_id"] is not None:
                    self._name = self.delta["seq_id"]
                else:
                    pass
            elif self.gamma is not None:
                if self.gamma.id is not None:
                    self._name = self.gamma.id
                elif self.gamma["sequence_id"] is not None:
                    self._name = self.gamma["sequence_id"]
                elif self.gamma["seq_id"] is not None:
                    self._name = self.gamma["seq_id"]
                else:
                    pass
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    # @property
    # def sample(self):
    #     if self._sample is None:
    #         slist = []
    #         if self.experiment is not None:
    #             slist.append(str(self.experiment))
    #         if self.group is not None:
    #             slist.append(str(self.group))
    #         if self.subject is not None:
    #             slist.append(str(self.subject))
    #         if self.timepoint is not None:
    #             slist.append(str(self.timepoint))
    #         if slist:
    #             self._sample = '|'.join(slist)
    #     return self._sample

    # @property
    # def subject(self):
    #     if self._subject is None:
    #         if self.heavy is not None and 'subject' in self.heavy.keys():
    #             self._subject = self.heavy['subject']
    #         elif self.light is not None and 'subject' in self.light.keys():
    #             self._subject = self.light['subject']
    #     return self._subject

    # @subject.setter
    # def subject(self, subject):
    #     self._subject = subject

    # @property
    # def group(self):
    #     if self._group is None:
    #         if self.heavy is not None and 'group' in self.heavy.keys():
    #             self._group = self.heavy['group']
    #         elif self.light is not None and 'group' in self.light.keys():
    #             self._group = self.light['group']
    #     return self._group

    # @group.setter
    # def group(self, group):
    #     self._group = group

    # @property
    # def experiment(self):
    #     if self._experiment is None:
    #         if self.heavy is not None and 'experiment' in self.heavy.keys():
    #             self._experiment = self.heavy['experiment']
    #         elif self.light is not None and 'experiment' in self.light.keys():
    #             self._experiment = self.light['experiment']
    #     return self._experiment

    # @experiment.setter
    # def experiment(self, experiment):
    #     self._experiment = experiment

    # @property
    # def timepoint(self):
    #     if self._timepoint is None:
    #         if self.heavy is not None and 'timepoint' in self.heavy.keys():
    #             self._timepoint = self.heavy['timepoint']
    #         elif self.light is not None and 'timepoint' in self.light.keys():
    #             self._timepoint = self.light['timepoint']
    #     return self._timepoint

    # @timepoint.setter
    # def timepoint(self, timepoint):
    #     self._timepoint = timepoint

    # def refine(self, heavy=True, light=True, species='human'):
    #     for seq in [s for s in [self.heavy, self.light] if s is not None]:
    #         try:
    #             self.remove_ambigs(seq)
    #             self._refine_v(seq, species)
    #             self._refine_j(seq, species)
    #             self._retranslate(seq)
    #         except:
    #             print('REFINEMENT FAILED: {}, {} chain'.format(seq['seq_id'], seq['chain']))
    #             print(traceback.format_exception_only(*sys.exc_info()[:2]))

    # @staticmethod
    # def remove_ambigs(seq):
    #     # fix Ns in the nucleotide sequence
    #     vdj = ''
    #     for s, g in zip(seq['vdj_nt'], seq['vdj_germ_nt']):
    #         if s.upper() == 'N':
    #             vdj += g
    #         else:
    #             vdj += s
    #     seq['vdj_nt'] = vdj
    #     # fix Xs in the amino acid sequence
    #     vdj = ''
    #     for s, g in zip(seq['vdj_aa'], seq['vdj_germ_aa']):
    #         if s.upper() == 'X':
    #             vdj += g
    #         else:
    #             vdj += s
    #     seq['vdj_aa'] = vdj

    # @staticmethod
    # def _refine_v(seq, species):
    #     '''
    #     Completes the 5' end of a a truncated sequence with germline nucleotides.
    #     Input is a MongoDB dict (seq) and the species.
    #     '''
    #     vgerm = germlines.get_germline(seq['v_gene']['full'], species)
    #     aln = global_alignment(seq['vdj_nt'], vgerm)
    #     prepend = ''
    #     for s, g in zip(aln.aligned_query, aln.aligned_target):
    #         if s != '-':
    #             break
    #         else:
    #             prepend += g
    #     seq['vdj_nt'] = prepend + seq['vdj_nt']

    # @staticmethod
    # def _refine_j(seq, species):
    #     '''
    #     Completes the 3' end of a a truncated sequence with germline nucleotides.
    #     Input is a MongoDB dict (seq) and the species.
    #     '''
    #     jgerm = germlines.get_germline(seq['j_gene']['full'], species)
    #     aln = global_alignment(seq['vdj_nt'], jgerm)
    #     append = ''
    #     for s, g in zip(aln.aligned_query[::-1], aln.aligned_target[::-1]):
    #         if s != '-':
    #             break
    #         else:
    #             append += g
    #     seq['vdj_nt'] = seq['vdj_nt'] + append[::-1]

    # @staticmethod
    # def _retranslate(seq):
    #     '''
    #     Retranslates a nucleotide sequence following refinement.
    #     Input is a Pair sequence (basically a dict of MongoDB output).
    #     '''
    #     if len(seq['vdj_nt']) % 3 != 0:
    #         trunc = len(seq['vdj_nt']) % 3
    #         seq['vdj_nt'] = seq['vdj_nt'][:-trunc]
    #     seq['vdj_aa'] = Seq(seq['vdj_nt']).translate()

    def fasta(
        self, name_field="sequence_id", sequence_field="sequence", append_chain=True
    ):
        """
        Returns the sequence pair as a fasta string. If the Pair object contains
        both heavy and light chain sequences, both will be returned as a single string.

        By default, the fasta string contains the 'vdj_nt' sequence for each chain. To change,
        use the <key> option to select an alternate sequence.

        By default, the chain (heavy or light) will be appended to the sequence name:

        >MySequence_heavy

        To just use the pair name (which will result in duplicate sequence names for Pair objects
        with both heavy and light chains), set <append_chain> to False.
        """
        fastas = []
        for s, chain in [(self.heavy, "heavy"), (self.light, "light")]:
            if s is not None:
                c = "_{}".format(chain) if append_chain else ""
                fastas.append(">{}{}\n{}".format(s[name_field], c, s[sequence_field]))
        return "\n".join(fastas)

    def summarize(self, annotation_format="airr"):
        try:
            if annotation_format == "airr":
                v_key = "v_gene"
                d_key = "d_gene"
                j_key = "j_gene"
                junc_key = "junction_aa"
                vident_key = "v_identity"
                isotype_key = "isotype"
            elif annotation_format == "json":
                v_key = "v_gene.gene"
                d_key = "d_gene.gene"
                j_key = "j_gene.gene"
                junc_key = "junc_aa"
                vident_key = "nt_identity.v"
                isotype_key = "isotype"
            vline = ""
            dline = ""
            jline = ""
            juncline = ""
            identline = ""
            isotypeline = ""
            if self.heavy is not None:
                vline += nested_dict_lookup(self.heavy, v_key.split("."))
                dline += nested_dict_lookup(self.heavy, d_key.split("."))
                jline += nested_dict_lookup(self.heavy, j_key.split("."))
                juncline += nested_dict_lookup(self.heavy, junc_key.split("."))
                identline += nested_dict_lookup(self.heavy, vident_key.split("."))
                isotypeline += nested_dict_lookup(self.heavy, isotype_key.split("."))
            pad = (
                max(
                    [
                        len(l)
                        for l in [vline, dline, jline, juncline, identline, isotypeline]
                    ]
                )
                + 4
            )
            vline += " " * (pad - len(vline))
            dline += " " * (pad - len(dline))
            jline += " " * (pad - len(jline))
            juncline += " " * (pad - len(juncline))
            identline += " " * (pad - len(identline))
            isotypeline += " " * (pad - len(isotypeline))
            if self.light is not None:
                vline += nested_dict_lookup(self.light, v_key.split("."))
                jline += nested_dict_lookup(self.light, j_key.split("."))
                juncline += nested_dict_lookup(self.light, junc_key.split("."))
                identline += nested_dict_lookup(self.light, vident_key.split("."))

            header_len = max(
                [
                    len(l)
                    for l in [vline, dline, jline, juncline, identline, isotypeline]
                ]
            )
            header_spaces = int((header_len - len(self.name)) / 2)
            print(f"{' ' * header_spaces}{self.name}")
            print("-" * header_len)
            for l in [vline, dline, jline, juncline, identline, isotypeline]:
                print(l)
            print("")
        except:
            return

    def _chain_selector(self, seqs):
        if len(seqs) == 0:
            return None
        if all(["umis" in s for s in seqs]):
            sorted_seqs = sorted(seqs, key=lambda x: x["umis"], reverse=True)
            return sorted_seqs[0]
        else:
            return seqs[0]


def get_pairs(
    db,
    collection,
    experiment=None,
    subject=None,
    group=None,
    name="seq_id",
    delim=None,
    delim_occurance=1,
    pairs_only=False,
    chain_selection_func=None,
):
    """
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
    """
    match = {}
    if subject is not None:
        if type(subject) in (list, tuple):
            match["subject"] = {"$in": subject}
        elif type(subject) in STR_TYPES:
            match["subject"] = subject
    if group is not None:
        if type(group) in (list, tuple):
            match["group"] = {"$in": group}
        elif type(group) in STR_TYPES:
            match["group"] = group
    if experiment is not None:
        if type(experiment) in (list, tuple):
            match["experiment"] = {"$in": experiment}
        elif type(experiment) in STR_TYPES:
            match["experiment"] = experiment
    seqs = list(db[collection].find(match))
    return assign_pairs(
        seqs,
        name=name,
        delim=delim,
        delim_occurance=delim_occurance,
        pairs_only=pairs_only,
        chain_selection_func=chain_selection_func,
    )


def assign_pairs(
    seqs,
    id_key="sequence_id",
    delim=None,
    delim_occurance=1,
    pairs_only=False,
    chain_selection_func=None,
    tenx_annot_file=None,
):
    """
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
    """
    # add UMIs to the sequence annotations if a 10xG annotations file is provided
    if tenx_annot_file is not None:
        annots = {}
        umis = {}
        with open(tenx_annot_file) as f:
            reader = csv.DictReader(f)
            for r in reader:
                key = "consensus_id" if "consensus_id" in r else "contig_id"
                umis[r[key]] = r["umis"]
                annots[r[key]] = r
        for s in seqs:
            s["umis"] = int(umis.get(s[id_key], 0))
            s["tenx_annots"] = annots.get(s[id_key], {})

    # identify pairs
    pdict = {}
    for s in seqs:
        if delim is not None:
            pname = delim.join(s[id_key].split(delim)[:delim_occurance])
        else:
            pname = s[id_key]
        if pname not in pdict:
            pdict[pname] = [
                s,
            ]
        else:
            pdict[pname].append(s)
    pairs = [
        Pair(p, name=n, chain_selection_func=chain_selection_func)
        for n, p in pdict.items()
    ]
    if pairs_only:
        pairs = [p for p in pairs if p.is_pair]
    return pairs


# def deduplicate(pairs, aa=False, ignore_primer_regions=False):
#     '''
#     Removes duplicate sequences from a list of Pair objects.

#     If a Pair has heavy and light chains, both chains must identically match heavy and light chains
#     from another Pair to be considered a duplicate. If a Pair has only a single chain,
#     identical matches to that chain will cause the single chain Pair to be considered a duplicate,
#     even if the comparison Pair has both chains.

#     Note that identical sequences are identified by simple string comparison, so sequences of
#     different length that are identical over the entirety of the shorter sequence are not
#     considered duplicates.

#     By default, comparison is made on the nucleotide sequence. To use the amino acid sequence instead,
#     set aa=True.
#     '''
#     nr_pairs = []
#     just_pairs = [p for p in pairs if p.is_pair]
#     single_chains = [p for p in pairs if not p.is_pair]
#     _pairs = just_pairs + single_chains
#     for p in _pairs:
#         duplicates = []
#         for nr in nr_pairs:
#             identical = True
#             vdj = 'vdj_aa' if aa else 'vdj_nt'
#             offset = 4 if aa else 12
#             if p.heavy is not None:
#                 if nr.heavy is None:
#                     identical = False
#                 else:
#                     heavy = p.heavy[vdj][offset:-offset] if ignore_primer_regions else p.heavy[vdj]
#                     nr_heavy = nr.heavy[vdj][offset:-offset] if ignore_primer_regions else nr.heavy[vdj]
#                     if heavy != nr_heavy:
#                         identical = False
#             if p.light is not None:
#                 if nr.light is None:
#                     identical = False
#                 else:
#                     light = p.light[vdj][offset:-offset] if ignore_primer_regions else p.light[vdj]
#                     nr_light = nr.light[vdj][offset:-offset] if ignore_primer_regions else nr.light[vdj]
#                     if light != nr_light:
#                         identical = False
#             duplicates.append(identical)
#         if any(duplicates):
#             continue
#         else:
#             nr_pairs.append(p)
#     return nr_pairs


# def refine(pairs, heavy=True, light=True, species='human'):
#     refined_pairs = copy.deepcopy(pairs)
#     for p in refined_pairs:
#         p.refine(heavy, light, species)
#     return refined_pairs


def vrc01_like(pair, loose=False):
    if loose:
        return loose_vrc01_like(pair)
    else:
        return strict_vrc01_like(pair)


def strict_vrc01_like(pair):
    if any([pair.heavy is None, pair.light is None]):
        vrc01_like = False
    else:
        try:
            # abstar's JSON output format
            vrc01_like = all(
                [pair.heavy["v_gene"]["gene"] == "IGHV1-2", pair.light["cdr3_len"] == 5]
            )
        except KeyError:
            try:
                # AIRR format
                vrc01_like = all(
                    [
                        pair.heavy["v_call"] == "IGHV1-2",
                        pair.light["junction_aa_length"] == 7,
                    ]
                )
            except KeyError:
                # unknown format
                vrc01_like = None
    return vrc01_like


def loose_vrc01_like(pair):
    try:
        # abstar's JSON output format
        if any([l["cdr3_len"] == 5 for l in pair.lights]) and any(
            [h["v_gene"]["gene"] == "IGHV1-2" for h in pair.heavies]
        ):
            loose_vrc01_like = True
        else:
            loose_vrc01_like = False
    except KeyError:
        try:
            # AIRR format
            if any([l["junction_aa_length"] == 7 for l in pair.lights]) and any(
                [h["v_call"] == "IGHV1-2" for h in pair.heavies]
            ):
                loose_vrc01_like = True
            else:
                loose_vrc01_like = False
        except:
            # unknown format
            loose_vrc01_like = None
    return loose_vrc01_like
