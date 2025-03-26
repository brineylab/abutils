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

import csv
import os
from typing import Callable, Iterable, Optional, Union

import pandas as pd
import polars as pl

from ..core.sequence import Sequence

# from ..utils.utilities import nested_dict_lookup


class Pair(object):
    """
    Holds a one or more ``Sequence`` objects, corresponding to a paired mAb or TCR.

    """

    def __init__(
        self,
        sequences: Iterable[Union[dict, Sequence]],
        name: Optional[str] = None,
        chain_selection_func: Optional[Callable] = None,
        properties: Optional[dict] = None,
    ):
        """
        Initialize a Pair object.

        Parameters
        ----------
        sequences : Iterable[Union[dict, Sequence]]
            A list of sequence objects, each containing sequence information.

        name : Optional[str], default=None
            The name of the pair.

        chain_selection_func : Optional[Callable], default=None
            A function that takes a list of sequences and orders them to determine
            the "correct" heavy and light chains in cases for which multiple heavy
            or light chains exist. If not provided, chains are prioritized in the order
            provided.

        properties : Optional[dict], default=None
            A dictionary of additional properties to add to the Pair object.

        """
        self._sequences = sequences
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
        if properties is not None:
            for k, v in properties.items():
                setattr(self, k, v)

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
            # check for the "locus" annotation
            if all([s["locus"] is not None for s in self._sequences]):
                if all(
                    [
                        s["locus"].upper() in ["IGH", "IGK", "IGL"]
                        for s in self._sequences
                    ]
                ):
                    self._receptor = "bcr"
                elif all(
                    [
                        s["locus"].upper() in ["TRA", "TRB", "TRD", "TRG"]
                        for s in self._sequences
                    ]
                ):
                    self._receptor = "tcr"
                else:
                    self._receptor = "unknown"

            # if locus isn't present, infer from V-gene assignment
            elif all([s["v_call"] is not None for s in self._sequences]):
                if all(
                    [
                        s["v_call"][:3].upper() in ["IGH", "IGK", "IGL"]
                        for s in self._sequences
                    ]
                ):
                    self._receptor = "bcr"
                elif all(
                    [
                        s["v_call"][:3].upper() in ["TRA", "TRB", "TRD", "TRG"]
                        for s in self._sequences
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
            self._heavies = [
                s for s in self._sequences if self._is_locus_type(s, ["IGH", "heavy"])
            ]
            # if all([s["chain"] is not None for s in self._seqs]):
            #     self._heavies = [s for s in self._seqs if s["chain"] == "heavy"]
            # elif all([s["locus"] is not None for s in self._seqs]):
            #     self._heavies = [s for s in self._seqs if s["locus"] == "IGH"]
            # elif all([s["v_gene"] is not None for s in self._seqs]):
            #     self._heavies = [
            #         s for s in self._seqs if s["v_gene"][:3].upper() == "IGH"
            #     ]
        return self._heavies

    @property
    def lights(self):
        if self._lights is None:
            self._lights = [
                s
                for s in self._sequences
                if self._is_locus_type(s, ["IGK", "IGL", "kappa", "lambda"])
            ]
            # if all([s["chain"] is not None for s in self._seqs]):
            #     self._lights = [
            #         s for s in self._seqs if s["chain"] in ["kappa", "lambda"]
            #     ]
            # elif all([s["locus"] is not None for s in self._seqs]):
            #     self._lights = [s for s in self._seqs if s["locus"] in ["IGK", "IGL"]]
            # elif all([s["v_gene"] is not None for s in self._seqs]):
            #     self._lights = [
            #         s for s in self._seqs if s["v_gene"][:3].upper() in ["IGK", "IGL"]
            #     ]
        return self._lights

    @property
    def alphas(self):
        if self._alphas is None:
            self._alphas = [
                s for s in self._sequences if self._is_locus_type(s, ["TRA", "alpha"])
            ]
            # if all([s["chain"] is not None for s in self._seqs]):
            #     self._alphas = [s for s in self._seqs if s["chain"] == "alpha"]
            # elif all([s["locus"] is not None for s in self._seqs]):
            #     self._alphas = [s for s in self._seqs if s["locus"] == "TRA"]
            # elif all([s["v_gene"] is not None for s in self._seqs]):
            #     self._alphas = [
            #         s for s in self._seqs if s["v_gene"][:3].upper() == "TRA"
            #     ]
        return self._alphas

    @property
    def betas(self):
        if self._betas is None:
            self._betas = [
                s for s in self._sequences if self._is_locus_type(s, ["TRB", "beta"])
            ]
            # if all([s["chain"] is not None for s in self._seqs]):
            #     self._betas = [s for s in self._seqs if s["chain"] == "beta"]
            # elif all([s["locus"] is not None for s in self._seqs]):
            #     self._betas = [s for s in self._seqs if s["locus"] == "TRB"]
            # elif all([s["v_gene"] is not None for s in self._seqs]):
            #     self._betas = [
            #         s for s in self._seqs if s["v_gene"][:3].upper() == "TRB"
            #     ]
        return self._betas

    @property
    def deltas(self):
        if self._deltas is None:
            self._deltas = [
                s for s in self._sequences if self._is_locus_type(s, ["TRD", "delta"])
            ]
            # if all([s["chain"] is not None for s in self._seqs]):
            #     self._deltas = [s for s in self._seqs if s["chain"] == "delta"]
            # elif all([s["locus"] is not None for s in self._seqs]):
            #     self._deltas = [s for s in self._seqs if s["locus"] == "TRD"]
            # elif all([s["v_gene"] is not None for s in self._seqs]):
            #     self._deltas = [
            #       s for s in self._seqs if s["v_gene"][:3].upper() == "TRD"
            #     ]
        return self._deltas

    @property
    def gammas(self):
        if self._gammas is None:
            self._gammas = [
                s for s in self._sequences if self._is_locus_type(s, ["TRG", "gamma"])
            ]
            # if all([s["chain"] is not None for s in self._seqs]):
            #     self._gammas = [s for s in self._seqs if s["chain"] == "gamma"]
            # elif all([s["locus"] is not None for s in self._seqs]):
            #     self._gammas = [s for s in self._seqs if s["locus"] == "TRG"]
            # elif all([s["v_gene"] is not None for s in self._seqs]):
            #     self._gammas = [
            #         s for s in self._seqs if s["v_gene"][:3].upper() == "TRG"
            #     ]
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

    @lineage.setter
    def lineage(self, lineage):
        self._lineage = lineage

    @property
    def name(self):
        if self._name is None:
            if self.heavy is not None:
                if self.heavy.id is not None:
                    self._name = self.heavy.id
                elif self.heavy["sequence_id"] is not None:
                    self._name = self.heavy["sequence_id"]
                # elif self.heavy["seq_id"] is not None:
                #     self._name = self.heavy["seq_id"]
                else:
                    pass
            elif self.light is not None:
                if self.light.id is not None:
                    self._name = self.light.id
                elif self.light["sequence_id"] is not None:
                    self._name = self.light["sequence_id"]
                # elif self.light["seq_id"] is not None:
                #     self._name = self.light["seq_id"]
                else:
                    pass
            elif self.alpha is not None:
                if self.alpha.id is not None:
                    self._name = self.alpha.id
                elif self.alpha["sequence_id"] is not None:
                    self._name = self.alpha["sequence_id"]
                # elif self.alpha["seq_id"] is not None:
                #     self._name = self.alpha["seq_id"]
                else:
                    pass
            elif self.beta is not None:
                if self.beta.id is not None:
                    self._name = self.beta.id
                elif self.beta["sequence_id"] is not None:
                    self._name = self.beta["sequence_id"]
                # elif self.beta["seq_id"] is not None:
                #     self._name = self.beta["seq_id"]
                else:
                    pass
            elif self.delta is not None:
                if self.delta.id is not None:
                    self._name = self.delta.id
                elif self.delta["sequence_id"] is not None:
                    self._name = self.delta["sequence_id"]
                # elif self.delta["seq_id"] is not None:
                #     self._name = self.delta["seq_id"]
                else:
                    pass
            elif self.gamma is not None:
                if self.gamma.id is not None:
                    self._name = self.gamma.id
                elif self.gamma["sequence_id"] is not None:
                    self._name = self.gamma["sequence_id"]
                # elif self.gamma["seq_id"] is not None:
                #     self._name = self.gamma["seq_id"]
                else:
                    pass
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    def _is_locus_type(self, seq: Sequence, locus_names: Iterable[str]) -> bool:
        """
        Check if the sequence is of the specified locus type.

        Parameters
        ----------
        seq : Sequence
            The sequence to check.

        locus_names : list[str]
            The list of locus names to check against.

        Returns
        -------
        bool
            True if the sequence is of the specified locus type, False otherwise.
        """
        if seq["locus"] in locus_names:
            return True
        elif seq["v_call"] is not None:
            if seq["v_call"][:3] in locus_names:
                return True
        return False

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
        self,
        name_field: str = "sequence_id",
        sequence_field: str = "sequence",
        append_chain: bool = True,
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

    def summarize(self, annotation_format: Optional[str] = None) -> None:
        """
        Prints a summary of the pair to the console.

        Parameters
        ----------
        annotation_format : str, optional
            No longer used, retained for legacy compatibility.

        Returns
        -------
        None
        """
        try:
            vline = ""
            dline = ""
            jline = ""
            juncline = ""
            identline = ""
            isotypeline = ""
            if self.heavy is not None:
                vline += self.heavy["v_gene"]
                dline += self.heavy["d_gene"]
                jline += self.heavy["j_gene"]
                juncline += self.heavy["junction_aa"]
                identline += self.heavy["v_identity"]
                isotypeline += self.heavy["isotype"]
            pad = (
                max(
                    [
                        len(line)
                        for line in [
                            vline,
                            dline,
                            jline,
                            juncline,
                            identline,
                            isotypeline,
                        ]
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
                vline += self.light["v_gene"]
                jline += self.light["j_gene"]
                juncline += self.light["junction_aa"]
                identline += self.light["v_identity"]
                isotypeline += self.light["isotype"]

            header_len = max(
                [
                    len(line)
                    for line in [vline, dline, jline, juncline, identline, isotypeline]
                ]
            )
            header_spaces = int((header_len - len(self.name)) / 2)
            print(f"{' ' * header_spaces}{self.name}")
            print("-" * header_len)
            for line in [vline, dline, jline, juncline, identline, isotypeline]:
                print(line)
            print("")
        except Exception:
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
        elif isinstance(subject, str):
            match["subject"] = subject
    if group is not None:
        if type(group) in (list, tuple):
            match["group"] = {"$in": group}
        elif isinstance(group, str):
            match["group"] = group
    if experiment is not None:
        if type(experiment) in (list, tuple):
            match["experiment"] = {"$in": experiment}
        elif isinstance(experiment, str):
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
    seqs: Iterable[Union[dict, Sequence]],
    id_key: str = "sequence_id",
    delim: Optional[str] = None,
    delim_occurance: int = 1,
    pairs_only: bool = False,
    chain_selection_func: Optional[Callable] = None,
    tenx_annot_file: Optional[str] = None,
) -> Iterable[Pair]:
    """
    Assigns sequences to the appropriate mAb pair, based on the sequence name.

    Parameters
    ----------
    seqs : Iterable[Union[dict, Sequence]]
        List of sequence objects, of the format returned by querying a MongoDB containing
        Abstar output.

    id_key : str, default="sequence_id"
        The dict key of the field to be used to group the sequences into pairs.
        Default is 'seq_id'

    delim : Optional[str], default=None
        An optional delimiter used to truncate the contents of the ::name:: field.
        Default is None, which results in no name truncation.

    delim_occurance : int, default=1
        The occurance of the delimiter at which to trim. Trimming is performed
        as delim.join(name.split(delim)[:delim_occurance]), so setting delim_occurance to -1 will
        trucate after the last occurance of delim. Default is 1.

    pairs_only : bool, default=False
        Setting to True results in only truly paired sequences (pair.is_pair == True)
        will be returned. Default is False.

    chain_selection_func : Optional[Callable], default=None
        A function that takes a list of sequences and returns a single sequence.
        Default is None, which results in the first sequence in the list being returned.

    tenx_annot_file : Optional[str], default=None
        A path to a 10x Genomics annotations file. If provided, the UMIs and 10x
        annotations will be added to the sequence annotations.

    Returns
    -------
    pairs : Iterable[Pair]
        A list of ``Pair`` objects, one for each mAb pair.

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


def pairs_to_polars(
    pairs: Iterable[Pair],
    annotations: Optional[Iterable[str]] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> pl.DataFrame:
    """
    Converts a list of ``Pair`` objects to a Polars DataFrame.

    Parameters
    ----------
    pairs : Iterable[Pair]
        List of ``Pair`` objects to be converted to a Polars DataFrame. Required.

    annotations : list, default=None
        A list of annotation fields to be included in the Polars DataFrame. The difference
        between ``annotations`` and ``columns`` is that ``annotations`` refers to fields
        in the heavy/light chain annotations (like ``"sequence"``, ``"v_gene"``, etc.),
        while ``columns`` refers to fields in the Polars DataFrame (like ``"sequence:0"``,
        ``"sequence:1"``, ``"name"``, etc.).

    columns : list, default=None
        A list of fields to be retained in the output Polars DataFrame. Fields must be column
        names in the input file. The difference between ``annotations`` and ``columns`` is
        that ``annotations`` refers to fields in the heavy/light chain annotations (like
        ``"sequence"``, ``"v_gene"``, etc.), while ``columns`` refers to fields in the
        Polars DataFrame (like ``"sequence:0"``, ``"sequence:1"``, ``"name"``, etc.).

    properties : list, default=None
        A list of properties to be included in the Polars DataFrame. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects.

    drop_na_columns : bool, default=True
        If ``True``, columns with all ``NaN`` values will be dropped from the Polars DataFrame.
        Default is ``True``.

    order : list, default=None
        A list of fields in the order they should appear in the Polars DataFrame.

    exclude : str or list, default=None
        Field or list of fields to be excluded from the Polars DataFrame.

    leading : str or list, default=None
        Field or list of fields to appear first in the Polars DataFrame. Supercedes ``order``, so
        if both are provided, fields in ``leading`` will appear first in the Polars DataFrame and
        remaining fields will appear in the order provided in ``order``.

    """
    data = []
    # read Pair data
    for p in pairs:
        d = {"name": p.name}
        # heavy chain
        if p.heavy is not None:
            if not p.heavy.annotations:
                d.update({"sequence_id:0": p.heavy.id, "sequence:0": p.heavy.sequence})
            else:
                if annotations is not None:
                    keys = [a for a in annotations if a in p.heavy.annotations]
                else:
                    keys = list(p.heavy.annotations.keys())
                d.update({f"{k}:0": p.heavy[k] for k in keys})
            if sequence_properties is not None:
                for prop in sequence_properties:
                    try:
                        d[f"{prop}:0"] = getattr(p.heavy, prop)
                    except AttributeError:
                        pass
        # light chain
        if p.light is not None:
            if not p.light.annotations:
                d.update({"sequence_id:1": p.light.id, "sequence:1": p.light.sequence})
            else:
                if annotations is not None:
                    keys = [a for a in annotations if a in p.light.annotations]
                else:
                    keys = list(p.light.annotations.keys())
                d.update({f"{k}:1": p.light[k] for k in keys})
            if sequence_properties is not None:
                for prop in sequence_properties:
                    try:
                        d[f"{prop}:1"] = getattr(p.light, prop)
                    except AttributeError:
                        pass
        # properties
        if properties is not None:
            for prop in properties:
                try:
                    d[prop] = getattr(p, prop)
                except AttributeError:
                    continue
        data.append(d)

    # build dataframe
    df = pl.DataFrame(data, infer_schema_length=len(data))

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

    if columns is not None:
        cols = [c for c in columns if c in df.columns]
        df = df.select(cols)

    return df


def pairs_from_polars(
    df: Union[pl.DataFrame, pl.LazyFrame],
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Pair]:
    """
    Reads a Polars DataFrame and returns ``Pair`` objects.

    Parameters
    ----------
    df : Union[pl.DataFrame, pl.LazyFrame]
        The input Polars DataFrame. Required.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'v_gene:0'`` field is not ``'IGHV1-2'``:

        .. code-block:: python

            {'v_gene:0': 'IGHV1-2'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Pair``
        objects. Fields must be column names in the input file.

    id_key : str, default="sequence_id"
        Name of the annotation field containing the sequence ID. Used to
        populate the ``Sequence.id`` property.

    sequence_key : str, default="sequence"
        Name of the annotation field containg the sequence. Used to
        populate the ``Sequence.sequence`` property.


    Returns
    -------
    pairs : list of ``Pairs``

    """
    if match is None:
        match = {}
    if fields is not None:
        if id_key not in fields:
            fields.append(id_key)
        if sequence_key not in fields:
            fields.append(sequence_key)
    pairs = []
    if isinstance(df, pl.LazyFrame):
        df = df.collect()
    for r in df.iter_rows(named=True):
        try:
            if all([r[k] == v for k, v in match.items()]):
                heavy = None
                light = None
                name = r.get("name", None)
                # heavy chain
                heavy_dict = {
                    k.split(":")[0]: v for k, v in r.items() if k.endswith(":0")
                }
                if fields is not None:
                    hc_fields = [f for f in fields if f in heavy_dict]
                    heavy_dict = {f: heavy_dict[f] for f in hc_fields}
                if any([v is not None for v in heavy_dict.values()]):
                    heavy = Sequence(heavy_dict, id_key=id_key, seq_key=sequence_key)
                # light chain
                light_dict = {
                    k.split(":")[0]: v for k, v in r.items() if k.endswith(":1")
                }
                if fields is not None:
                    lc_fields = [f for f in fields if f in light_dict]
                    light_dict = {f: light_dict[f] for f in lc_fields}
                if any([v is not None for v in light_dict.values()]):
                    light = Sequence(light_dict, id_key=id_key, seq_key=sequence_key)
                # extra fields:
                extra_fields = {}
                for k, v in r.items():
                    if (
                        not k.endswith(":0")
                        and not k.endswith(":1")
                        and k.lower() != "name"
                    ):
                        extra_fields[k] = v
                # pair
                sequences = [s for s in [heavy, light] if s is not None]
                pairs.append(
                    Pair(sequences=sequences, name=name, properties=extra_fields)
                )
        except KeyError:
            continue
    return pairs


def pairs_to_pandas(
    pairs: Iterable[Pair],
    annotations: Optional[Iterable[str]] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> pd.DataFrame:
    """
    Saves a list of ``Pair`` objects to a CSV file.

    Parameters
    ----------
    pairs : Iterable[Pair]
        List of ``Pair`` objects to be saved to a CSV file. Required.

    annotations : list, default=None
        A list of annotation fields to be included in the Polars DataFrame. The difference
        between ``annotations`` and ``columns`` is that ``annotations`` refers to fields
        in the heavy/light chain annotations (like ``"sequence"``, ``"v_gene"``, etc.),
        while ``columns`` refers to fields in the Polars DataFrame (like ``"sequence:0"``,
        ``"sequence:1"``, ``"name"``, etc.).

    columns : list, default=None
        A list of fields to be retained in the output Polars DataFrame. Fields must be column
        names in the input file. The difference between ``annotations`` and ``columns`` is
        that ``annotations`` refers to fields in the heavy/light chain annotations (like
        ``"sequence"``, ``"v_gene"``, etc.), while ``columns`` refers to fields in the
        Polars DataFrame (like ``"sequence:0"``, ``"sequence:1"``, ``"name"``, etc.).

    properties : list, default=None
        A list of properties to be included in the CSV file. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects.

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

    """
    df = pairs_to_polars(
        pairs,
        annotations=annotations,
        columns=columns,
        properties=properties,
        sequence_properties=sequence_properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )
    return df.to_pandas()


def pairs_from_pandas(
    df: pd.DataFrame,
    match: Optional[dict] = None,
    fields: Iterable = None,
    id_key: str = "sequence_id",
    sequence_key: str = "sequence",
) -> Iterable[Pair]:
    """
    Reads a Polars DataFrame and returns ``Pair`` objects.

    Parameters
    ----------
    df : pd.DataFrame
        The input pandas DataFrame. Required.

    match : dict, default=None
        A ``dict`` for filtering sequences from the input file.
        Sequences must match all conditions to be returned. For example,
        the following ``dict`` will filter out all sequences for which
        the ``'v_gene:0'`` field is not ``'IGHV1-2'``:

        .. code-block:: python

            {'v_gene:0': 'IGHV1-2'}

    fields : list, default=None
        A ``list`` of fields to be retained in the output ``Pair``
        objects. Fields must be column names in the input file.

    id_key : str, default="sequence_id"
        Name of the annotation field containing the sequence ID. Used to
        populate the ``Sequence.id`` property.

    sequence_key : str, default="sequence"
        Name of the annotation field containg the sequence. Used to
        populate the ``Sequence.sequence`` property.


    Returns
    -------
    pairs : list of ``Pairs``

    """
    polars_df = pl.from_pandas(df)
    return pairs_from_polars(
        df=polars_df,
        match=match,
        fields=fields,
        id_key=id_key,
        sequence_key=sequence_key,
    )


def pairs_to_csv(
    pairs: Iterable[Pair],
    csv_file: Optional[str] = None,
    separator: str = ",",
    header: bool = True,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Pair`` objects to a CSV file.

    Parameters
    ----------
    pairs : Iterable[Pair]
        List of ``Pair`` objects to be saved to a CSV file. Required.

    csv_file : str
        Path to the output CSV file. If not provided, a ``polars.DataFrame``
        containing the CSV data will be returned.

    separator : str, default=","
        Column separator. Default is ``","``.

    header : bool, default=True
        If ``True``, the CSV file will contain a header row. Default is ``True``.

    columns : list, default=None
        A list of fields to be retained in the output CSV file. Fields must be column
        names in the input file.

    properties : list, default=None
        A list of properties to be included in the CSV file. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects.

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

    """
    df = pairs_to_polars(
        pairs,
        columns=columns,
        properties=properties,
        sequence_properties=sequence_properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )
    df.write_csv(csv_file, separator=separator, include_header=header)


def pairs_to_parquet(
    pairs: Iterable[Pair],
    parquet_file: Optional[str] = None,
    columns: Optional[Iterable] = None,
    properties: Optional[Iterable[str]] = None,
    sequence_properties: Optional[Iterable[str]] = None,
    drop_na_columns: bool = True,
    order: Optional[Iterable[str]] = None,
    exclude: Optional[Union[str, Iterable[str]]] = None,
    leading: Optional[Union[str, Iterable[str]]] = None,
) -> None:
    """
    Saves a list of ``Pair`` objects to a Parquet file.

    Parameters
    ----------
    pairs : Iterable[Pair]
        List of ``Pair`` objects to be saved to a Parquet file. Required.

    parquet_file : str
        Path to the output Parquet file. If not provided, a ``polars.DataFrame``
        containing the Parquet data will be returned.

    columns : list, default=None
        A list of fields to be retained in the output Parquet file. Fields must be column
        names in the input file.

    properties : list, default=None
        A list of properties to be included in the Parquet file. If not provided, everything
        in the ``annotations`` field of each heavy/light chain will be included.

    sequence_properties : list, default=None
        A list of sequence properties to be included. Differs from ``properties``, which
        refers to properties of the ``Pair`` object. These properties are those of the
        heavy/light ``Sequence`` objects.

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

    """
    df = pairs_to_polars(
        pairs,
        columns=columns,
        properties=properties,
        sequence_properties=sequence_properties,
        drop_na_columns=drop_na_columns,
        order=order,
        exclude=exclude,
        leading=leading,
    )
    df.write_parquet(parquet_file)


# def pairs_to_fasta(
#     pairs: Iterable[Pair],
#     fasta_file: Optional[str] = None,
#     id_key: str = "sequence_id",
#     sequence_key: str = "sequence",
#     append_chain: bool = False,
# ) -> str:
#     """
#     Saves a list of ``Pair`` objects to a FASTA-formatted file.

#     Parameters
#     ----------
#     pairs : Iterable[Pair]
#         List of ``Pair`` objects to be saved to a FASTA-formatted file. Required.

#     fasta_file : str
#         Path to the output FASTA-formatted file. If not provided, the FASTA-formatted
#         sequences will be returned as a string.

#     Returns
#     -------
#     fasta_str : str
#         FASTA-formatted sequences as a string if ``fasta_file`` is not provided. If
#         ``fasta_file`` is provided, the path to the output file is returned.

#     """
#     fastas = [
#         p.to_fasta(id_key=id_key, sequence_key=sequence_key, append_chain=append_chain)
#         for p in pairs
#     ]
#     if fasta_file is None:
#         return "\n".join(fastas)
#     else:
#         from ..io import make_dir

#         make_dir(os.path.dirname(fasta_file))
#         with open(fasta_file, "w") as f:
#             f.write("\n".join(fastas))
#         return fasta_file


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
