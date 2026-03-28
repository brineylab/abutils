#!/usr/bin/env python
# filename: lineage.py

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


import math
import subprocess as sp
from collections import Counter
from collections.abc import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import prettytable as pt
from matplotlib.colors import ListedColormap
from natsort import natsorted

from ..io import from_polars, to_polars
from ..tools.alignment import muscle, semiglobal_alignment
from ..tools.cluster import cluster
from ..utils.decorators import lazy_property
from .pair import Pair
from .sequence import Sequence, translate


class Lineage:
    """
    Methods for manipulating an antibody lineage.


    INPUTS
    ------
    pairs: a list of one or more vaxtools.utils.pair.Pair objects


    PROPERTIES
    ----------
    name: the Clonify ID of the lineage (if Clonify was used to assign lineages).
        if ['clonify']['id'] does not exist in any of the heavy chains, name will
        be None.

    just_pairs: a list containing all lineage Pair objects that have both heavy
        and light chains.

    verified_pairs: a list containing all lineage Pair objects that contain a
        verified heavy/light pair.

    heavies: a list of all lineage Pair objects with a heavy chain (whether or not
        they also have a light chain)

    lights: a list of all lineage Pair objects with a light chain (whether or not
        they also have a heavy chain)

    uca: returns a Pair objecct of the unmutated common ancestor for the lineage.
        The uca is computed by taking the germline V(D)J regions, plus the N-addition
        region(s) from the least mutated heavy or light chain. If the lineage contains
        both heavy and light chains, the returned Pair will have ucas for both chains.

    subject: if any of the Pair objects contains a 'subject' property, this will be returned.
        If all 'subject' properties are the same, the return value will be a string. If there are
        multiple different 'subject' properties, the return value will be a list. If no Pairs
        have a 'subject' property, None will be returned.

    experiment: if any of the Pair objects contains a 'experiment' property, this will be returned.
        If all 'experiment' properties are the same, the return value will be a string. If there are
        multiple different 'experiment' properties, the return value will be a list. If no Pairs
        have a 'experiment' property, None will be returned.

    group: if any of the Pair objects contains a 'group' property, this will be returned.
        If all 'group' properties are the same, the return value will be a string. If there are
        multiple different 'group' properties, the return value will be a list. If no Pairs
        have a 'group' property, None will be returned.

    timepoint: if any of the Pair objects contains a 'timepoint' property, this will be returned.
        If all 'timepoint' properties are the same, the return value will be a string. If there are
        multiple different 'timepoint' properties, the return value will be a list. If no Pairs
        have a 'timepoint' property, None will be returned.


    METHODS
    -------

    Lineage.phylogeny()

        Inputs
        ------
        project_dir: directory for alignment, tree and figure files (required)

        aln_file: if an alignment file has already been computed, passing the file path
            will force phylogeny() to use this alignment instead of computing a new one.

        tree_file: if a tree file has already been computed, passing the file path
            will force phylogeny() to use this tree instead of computing a new one.

        aa: if True, build an alignment based on amino acid sequences. Default is False.

        root: provide a sequence (string) to root the tree. Default is to use the UCA.

        colors: a dict that maps sequences to colors (as hex values). If used without
            orders or order_function, keys should be Pair names and values should be hex
            strings for all Pairs to be colored. Any Pair not present in the dict will
            be colored black. If used with orders or order_function, the keys should be
            orders and the values should be hex strings. As before, any order values
            not present in the dict will be colored black.

    """

    def __init__(self, pairs: Iterable[Pair], name: str | None = None):
        """
        Initialize a ``Lineage`` object.

        Parameters
        ----------
        pairs : Iterable[Pair]
            An iterable of ``Pair`` objects

        name : Optional[str]
            The name of the lineage (optional)

        """
        self.pairs = pairs
        self.rmp_threshold = 0.9
        self.rmp_alt_allowed_mismatches = 1

    def __contains__(self, item):
        if item in self.pair_dict:
            return True
        return False

    def __getitem__(self, key):
        return self.pair_dict.get(key, None)

    def __setitem__(self, key, val):
        self.pair_dict[key] = val

    def __iter__(self):
        return iter(self.pairs)

    @lazy_property
    def pair_dict(self):
        return {p.name: p for p in self.pairs}

    @lazy_property
    def name(self):
        """
        Returns the lineage name, or ``None`` if the name cannot be found.
        """
        if any([hasattr(p, "lineage") for p in self.pairs]):
            lineage_names = [p.lineage for p in self.pairs if hasattr(p, "lineage")]
        else:
            lineage_names = [
                p.heavy["lineage"] for p in self.heavies if "lineage" in p.heavy
            ]
        if len(lineage_names) > 0:
            return lineage_names[0]
        return None

    @lazy_property
    def just_pairs(self):
        """
        Returns all lineage Pair objects that contain both
        heavy and light chains.
        """
        return [p for p in self.pairs if p.is_pair]

    @lazy_property
    def verified_pairs(self):
        """
        Returns all lineage Pair objects that contain verified pairings.
        """
        if not hasattr(self.just_pairs[0], "verified"):
            self.verify_light_chains()
        return [p for p in self.just_pairs if p.verified]

    @lazy_property
    def heavies(self):
        """
        Returns all lineage Pair objects that contain a heavy chain.
        """
        return [p for p in self.pairs if p.heavy is not None]

    @lazy_property
    def lights(self):
        """
        Returns all lineage Pair objects that contain a light chain.
        """
        return [p for p in self.pairs if p.light is not None]

    @lazy_property
    def uca(self):
        """
        Calculates the unmutated common ancestor (UCA) of the lineage.

        The UCA is computed by using the germline V(D)J genes as well
        as the N-addition regions of the least mutated sequence.

        If both heavy and light chains exist in the lineage, the UCA
        will be computed for both chains. If not, the UCA will only
        be computed on the chains that exist in the lineage.

        Returns
        -------
        A VaxTools Pair object containing the UCA.
        """
        return self._calculate_uca()

    @lazy_property
    def rmp(self):
        return self._calculate_recalled_memory_precursor()

    @lazy_property
    def subject(self):
        return self._get_metadata("subject")

    @lazy_property
    def group(self):
        return self._get_metadata("group")

    @lazy_property
    def experiment(self):
        return self._get_metadata("experiment")

    @lazy_property
    def timepoints(self):
        if self.heavies:
            _timepoints = [p.heavy._mongo.get("timepoint", None) for p in self.heavies]
            return list(set([g for g in _timepoints if g is not None]))
        if self.lights:
            _timepoints = [p.light._mongo.get("timepoint", None) for p in self.lights]
            return list(set([g for g in _timepoints if g is not None]))
        return []

    @lazy_property
    def has_insertion(self):
        ins = ["v_ins" in p.heavy for p in self.heavies]
        ins += ["v_ins" in p.light for p in self.lights]
        return True if any(ins) else False

    @lazy_property
    def has_deletion(self):
        dels = ["v_del" in p.heavy for p in self.heavies]
        dels += ["v_del" in p.light for p in self.lights]
        return True if any(dels) else False

    @lazy_property
    def has_indel(self):
        if any([self.has_insertion, self.has_deletion]):
            return True
        return False

    def size(self, pairs_only: bool = False) -> int:
        """
        Calculate the size of the lineage.

        Parameters
        ----------
        pairs_only: bool, default=False
            If True, only count pairs with both heavy and light chains.

        Returns
        -------
        int
            The size of the lineage.

        """
        if pairs_only:
            return len(self.just_pairs)
        else:
            return len(self.heavies)

    def verify_light_chains(self, threshold: float = 0.9):
        """
        Clusters the light chains to identify potentially spurious (non-lineage)
        pairings. Following clustering, all pairs in the largest light chain
        cluster are assumed to be correctly paired. For each of those pairs,
        the ``verified`` attribute is set to True. For pairs not in the largest
        light chain cluster, the ``verified`` attribute is set to False.

        Parameters
        ----------
        threshold: float, default=0.9
            The clustering threshold.

        """
        lseqs = [l.light for l in self.lights]
        clusters = cluster(lseqs, threshold=threshold)
        clusters.sort(key=lambda x: x.size, reverse=True)
        verified_ids = clusters[0].ids
        for p in self.lights:
            p.verified = True if p.name in verified_ids else False

    def dot_alignment(
        self,
        seq_field: str = "sequence",
        name_field: str | None = None,
        uca: Sequence | None = None,
        chain: str = "heavy",
        uca_name: str = "UCA",
        as_fasta: bool = False,
        just_alignment: bool = False,
    ):
        """
        Returns a multiple sequence alignment of all lineage sequence with the UCA
        where matches to the UCA are shown as dots and mismatches are shown as the
        mismatched residue.

        Parameters
        ----------
        seq_field: str, default='sequence'
            The sequence field to be used for alignment.

        name_field: Optional[str], default=None
            The field used for the sequence name. If None, the ``Pair`` name will be used.

        chain: str, default='heavy'
            Either 'heavy' or 'light'.

        uca: Optional[Sequence], default=None
            The UCA sequence. If not provided, the UCA will be computed.

        uca_name: str, default='UCA'
            The name of the UCA sequence.

        as_fasta: bool, default=False
            If True, return the alignment as a FASTA string.

        just_alignment: bool, default=False
            If True, return the alignment as a list of strings.

        Returns
        -------
        The dot alignment (string or list of strings)

        """
        if uca is None:
            uca = self.uca.heavy if chain == "heavy" else self.uca.light
        uca.id = "UCA"
        if chain == "heavy":
            sequences = [p.heavy for p in self.heavies if seq_field in p.heavy]
            if name_field != "seq_id":
                uca[name_field] = uca["seq_id"]
            sequences.append(uca)
            seqs = [(s[name_field], s[seq_field]) for s in sequences]
        else:
            sequences = [p.light for p in self.lights if seq_field in p.light]
            if name_field != "seq_id":
                uca[name_field] = uca["seq_id"]
            sequences.append(uca)
            seqs = [(s[name_field], s[seq_field]) for s in sequences]
        aln = muscle(seqs)
        g_aln = [a for a in aln if a.id == "UCA"][0]
        dots = [
            (uca_name, str(g_aln.seq)),
        ]
        for seq in [a for a in aln if a.id != "UCA"]:
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

    def pixel(self, seq_field="vdj_nt", figfile=None):
        if "aa" in seq_field:
            colors = _aa_pixel_colors
        else:
            colors = _nt_pixel_colors
        ckeys = sorted(colors.keys())
        res_vals = {r: v for r, v in zip(ckeys, range(len(ckeys)))}
        cmap_colors = [colors[res] for res in ckeys]
        cmap = ListedColormap(cmap_colors)
        data = []
        dot_alignment = self.dot_alignment(seq_field=seq_field, just_alignment=True)
        for dot in dot_alignment:
            data.append([res_vals.get(res.upper(), len(res_vals) + 1) for res in dot])
        mag = (int(math.log10(len(data[0]))) + int(math.log10(len(data)))) / 2
        x_dim = len(data[0]) / 10**mag
        y_dim = len(data) / 10**mag
        plt.figure(figsize=(x_dim, y_dim), dpi=100)
        plt.imshow(data, cmap=cmap, interpolation="none")
        plt.axis("off")
        if figfile is not None:
            plt.savefig(figfile, bbox_inches="tight", dpi=400)
            plt.close()
        else:
            plt.show()

    def _calculate_uca(self, paired_only=False):
        from abstar import run as run_abstar

        if paired_only:
            heavies = self.just_pairs
            lights = self.just_pairs
        else:
            heavies = self.heavies
            lights = self.lights
        # heavy chain UCA
        if len(heavies) >= 1:
            lmhc = sorted(
                heavies, key=lambda x: x.heavy["nt_identity"]["v"], reverse=True
            )[0].heavy
            hc_uca = run_abstar(("UCA", lmhc["vdj_germ_nt"]), isotype=False)
        else:
            hc_uca = None
        # light chain UCA
        if len(lights) >= 1:
            lmlc = sorted(
                lights, key=lambda x: x.light["nt_identity"]["v"], reverse=True
            )[0].light
            lc_uca = run_abstar(("UCA", lmlc["vdj_germ_nt"]), isotype=False)
        else:
            lc_uca = None
        return Pair([uca for uca in [hc_uca, lc_uca] if uca is not None])

    def _calculate_recalled_memory_precursor(self, paired_only=False):
        if paired_only:
            heavies = [p.heavy for p in self.just_pairs]
            lights = [p.light for p in self.just_pairs]
        else:
            heavies = [p.heavy for p in self.heavies]
            lights = [p.light for p in self.lights]
        if len(heavies) >= 1:
            hc_rmp = self._rmp(heavies + [self.uca.heavy])
        else:
            hc_rmp = None
        if len(lights) >= 1:
            lc_rmp = self._rmp(lights + [self.uca.light])
        else:
            lc_rmp = None
        return Pair([hc_rmp, lc_rmp])

    def _rmp(self, sequences):
        from abstar import run as run_abstar

        rmp = ""
        seqs = [(s["seq_id"], s["vdj_nt"]) for s in sequences]
        aln = muscle(seqs)
        g_aln = [a for a in aln if a.id == "UCA"][0]
        query_seqs = [str(a.seq) for a in aln if a.id != "UCA"]
        for i, g in enumerate(g_aln):
            qcounts = Counter([q[i] for q in query_seqs])
            qmax = sorted(qcounts.keys(), key=lambda x: qcounts[x], reverse=True)[0]
            qmax_fraction = float(qcounts[qmax]) / sum(qcounts.values())
            qmax_alt_mismatches = sum(qcounts.values()) - qcounts[qmax]
            if any(
                [
                    qmax_fraction >= self.rmp_threshold,
                    qmax_alt_mismatches <= self.rmp_alt_allowed_mismatches,
                ]
            ):
                rmp += qmax
            else:
                rmp += g
        return run_abstar(Sequence(rmp.replace("-", ""), id="RMP"))

    def _germline_field_map(self, field):
        fmap = {
            "fr1_nt": "fr1_germ_nt",
            "fr2_nt": "fr2_germ_nt",
            "fr3_nt": "fr3_germ_nt",
            "fr4_nt": "fr4_germ_nt",
            "cdr1_nt": "cdr1_germ_nt",
            "cdr2_nt": "cdr2_germ_nt",
            "fr1_aa": "fr1_germ_aa",
            "fr2_aa": "fr2_germ_aa",
            "fr3_aa": "fr3_germ_aa",
            "fr4_aa": "fr4_germ_aa",
            "cdr1_aa": "cdr1_germ_aa",
            "cdr2_aa": "cdr2_germ_aa",
            "vdj_nt": "vdj_germ_nt",
            "vdj_aa": "vdj_germ_aa",
        }
        return fmap.get(field.lower(), None)


class LineageSummary:
    """
    docstring for LineageSummary
    """

    def __init__(
        self,
        lineage,
        name: str | None = None,
        dot_alignment: bool = True,
        in_color: bool = True,
        pairs_only: bool = False,
        padding_width: int = 1,
    ):
        if isinstance(lineage, (pl.DataFrame, pl.LazyFrame)):
            if isinstance(lineage, pl.LazyFrame):
                lineage = lineage.collect()
            self.df = lineage
            self.pairs = from_polars(lineage)
        elif isinstance(lineage, pd.DataFrame):
            lineage = pl.from_pandas(lineage)
            self.df = lineage
            self.pairs = from_polars(lineage)
        elif all([isinstance(l, Pair) for l in lineage]):
            self.pairs = lineage
            self.df = to_polars(lineage)
        else:
            raise ValueError("invalid input type")
        self.lineage = lineage
        self.name = name

        self.heavies = [p.heavy for p in self.pairs if p.heavy is not None]
        self.lights = [p.light for p in self.pairs if p.light is not None]

        self.adata = None
        self.obs = None
        self.bcrs = None
        self.agbcs = None
        self.specificities = None
        self.extra_blocks = []
        self.dot_alignment = dot_alignment
        self.in_color = in_color
        self.pairs_only = pairs_only
        self.padding_width = padding_width

        self.COLOR_PREFIX = {
            "yellow": "\033[93m",
            "green": "\033[92m",
            "red": "\033[91m",
            "blue": "\033[94m",
            "pink": "\033[95m",
            "black": "\033[90m",
        }

        self.REGION_COLORS = {
            "cdr1": "red",
            "cdr2": "yellow",
            "junction_v": "green",
            "junction_j": "blue",
            "junction": "black",
        }

    def __repr__(self):
        return self._assemble_summary()

    def __str__(self):
        s = self._assemble_summary()
        return s.get_string()

    def show(self, in_color=None):
        s = self._assemble_summary(in_color=in_color)
        print(s)

    def to_string(self, in_color=False):
        s = self._assemble_summary(in_color=in_color)
        return s.get_string()

    def _assemble_summary(self, in_color=None):
        in_color = in_color if in_color is not None else self.in_color

        # group identical AA sequences
        identical_dict = {}
        for pair in self.pairs:
            aa = ""
            if pair.heavy is not None:
                aa += pair.heavy["sequence_aa"]
            if pair.light is not None:
                aa += pair.light["sequence_aa"]
            if aa not in identical_dict:
                identical_dict[aa] = []
            identical_dict[aa].append(pair)
        identical_groups = sorted(
            identical_dict.values(), key=lambda x: len(x), reverse=True
        )
        h_standards = [ig[0].heavy for ig in identical_groups]
        l_standards = [ig[0].light for ig in identical_groups]

        # initialize table
        x = pt.PrettyTable()

        # frequencies block
        freq_block = self._frequencies_block(identical_groups)
        x.add_column("frequencies", ["", freq_block], align="r")

        # # specificities_block
        # if self.specificities is not None:
        #     spec_block = self._specificities_block(identical_groups)
        #     x.add_column("specificities", ["", spec_block], align="r")

        # # AgBCs block
        # if self.agbcs is not None:
        #     agbc_block = self._agbc_block(identical_groups)
        #     x.add_column("agbcs", ["", agbc_block], align="r")

        # # extra blocks
        # for xblock in self.extra_blocks:
        #     extra_block = self._extra_block(xblock, identical_groups)
        #     x.add_column(xblock, ["", extra_block], align="r")

        # HC and LC blocks
        hgene_block = self._germline_block(self.heavies)
        hcdr_block = self._cdr_alignment_block(h_standards, color=in_color)
        x.add_column("heavy", [hgene_block, hcdr_block], align="l")

        lgene_block = self._germline_block(self.lights)
        lcdr_block = self._cdr_alignment_block(l_standards, color=in_color)
        x.add_column("light", [lgene_block, lcdr_block], align="l")

        # style
        x.set_style(pt.DEFAULT)
        x.header = False
        x.hrules = pt.ALL
        x.left_padding_width = self.padding_width
        x.right_padding_width = self.padding_width
        x.title = self.name
        return x

    def _frequencies_block(self, identical_groups):
        # id column
        ids = ["", "all"] + [f"{i}" for i in range(len(identical_groups))]

        # sequence count
        n_seqs = ["n", len(self.pairs)] + [len(ig) for ig in identical_groups]

        # isotypes
        isotypes = natsorted(set([h["c_gene"] for h in self.heavies]))
        isotype_dict = {}
        for i in isotypes:
            isotype_dict[i] = []
            for ig in identical_groups:
                heavies = [b.heavy for b in ig if b is not None]
                isotype_dict[i].append(sum([h["c_gene"] == i for h in heavies]))
            isotype_dict[i] = [i, sum(isotype_dict[i])] + isotype_dict[i]

        # build the block
        x = pt.PrettyTable()
        x.add_column("ids", ids)
        x.add_column("n_seqs", n_seqs)
        for i in isotypes:
            x.add_column(i, isotype_dict[i])
        x.set_style(pt.PLAIN_COLUMNS)
        x.header = False
        x.left_padding_width = 0
        x.right_padding_width = self.padding_width
        x.align = "r"
        table_string = x.get_string()
        unpadded = "\n".join([l.strip() for l in table_string.split("\n")])
        return unpadded

    # def _specificities_block(self, identical_groups):
    #     specificities_dict = {}
    #     for s in self.specificities:
    #         counts = []
    #         for ig in identical_groups:
    #             names = [b.name for b in ig]
    #             positives = sum(self.obs[f"is_{s}"][names])
    #             counts.append(positives)
    #         specificities_dict[s] = [s, sum(counts)] + counts
    #     # build the block
    #     x = pt.PrettyTable()
    #     for s in self.specificities:
    #         x.add_column(s, specificities_dict[s])
    #     x.set_style(pt.PLAIN_COLUMNS)
    #     x.header = False
    #     x.left_padding_width = 0
    #     x.right_padding_width = self.padding_width
    #     x.align = "r"
    #     table_string = x.get_string()
    #     unpadded = "\n".join([l.strip() for l in table_string.split("\n")])
    #     return unpadded

    # def _agbc_block(self, identical_groups):
    #     agbc_dict = {}
    #     for a in self.agbcs:
    #         vals = []
    #         for ig in identical_groups:
    #             names = [b.name for b in ig]
    #             mean = np.mean(self.obs[f"{a}"][names])
    #             std = np.std(self.obs[f"{a}"][names])
    #             vals.append(f"{mean:.1f} ({std:.1f})")
    #         agbc_dict[a] = [a, ""] + vals
    #     # build the block
    #     x = pt.PrettyTable()
    #     for a in self.agbcs:
    #         x.add_column(a, agbc_dict[a])
    #     x.set_style(pt.PLAIN_COLUMNS)
    #     x.header = False
    #     x.left_padding_width = 0
    #     x.right_padding_width = self.padding_width
    #     x.align = "r"
    #     table_string = x.get_string()
    #     unpadded = "\n".join([l.strip() for l in table_string.split("\n")])
    #     return unpadded

    # def _extra_block(self, column, identical_groups):
    #     categories = natsorted(self.adata.obs[column].unique())
    #     extra_dict = {}
    #     for c in categories:
    #         extra_dict[c] = []
    #         for ig in identical_groups:
    #             names = [b.name for b in ig]
    #             _adata = self.adata[names]
    #             extra_dict[c].append(_adata[_adata.obs[column] == c].shape[0])
    #         extra_dict[c] = [c, sum(extra_dict[c])] + extra_dict[c]
    #     # build the block
    #     x = pt.PrettyTable()
    #     for c in categories:
    #         x.add_column(c, extra_dict[c])
    #     x.set_style(pt.PLAIN_COLUMNS)
    #     x.header = False
    #     x.left_padding_width = 0
    #     x.right_padding_width = self.padding_width
    #     x.align = "r"
    #     table_string = x.get_string()
    #     unpadded = "\n".join([l.strip() for l in table_string.split("\n")])
    #     return unpadded

    def _germline_block(self, sequences):
        block_strings = []

        # chain
        sequences = [s for s in sequences if s is not None]
        if all([s["locus"] == "IGH" for s in sequences]):
            chain = "HEAVY CHAIN\n"
        else:
            chain = "LIGHT CHAIN\n"
        # germline genes and CDR3 length
        v_call = Counter([s["v_call"] for s in sequences]).most_common(1)[0][0]
        if any([s["d_call"] for s in sequences]):
            d_call = Counter(
                [s["d_call"] for s in sequences if s["d_call"]]
            ).most_common(1)[0][0]
            dcount = len(list(set([s["d_call"] for s in sequences])))
            if len(list(set([s["d_call"] for s in sequences]))) > 1:
                d_call = f"{d_call}{'.' * (dcount-1)}"
        else:
            d_call = None
        j_call = Counter([s["j_call"] for s in sequences]).most_common(1)[0][0]
        cdr3_length = Counter([s["cdr3_length"] for s in sequences]).most_common(1)[0][
            0
        ]
        if d_call is not None:
            block_strings.append(
                f"{chain}{v_call} | {d_call} | {j_call} | {cdr3_length}"
            )
        else:
            block_strings.append(f"{chain}{v_call} | {j_call} | {cdr3_length}")

        # identities
        ident_string = "mutation: "
        nt_identities = [100.0 * (1 - float(s["v_identity"])) for s in sequences]
        min_nt = round(np.min(nt_identities), 1)
        max_nt = round(np.max(nt_identities), 1)
        if f"{min_nt:.0f}" == f"{max_nt:.0f}":
            ident_string += f"{min_nt:.0f}% nt"
        else:
            ident_string += f"{min_nt:.0f}-{max_nt:.0f}% nt"
        aa_identities = [100.0 * (1 - float(s["v_identity_aa"])) for s in sequences]
        mean_aa = round(np.mean(aa_identities), 1)
        min_aa = round(np.min(aa_identities), 1)
        max_aa = round(np.max(aa_identities), 1)
        ident_string += " | "
        if f"{min_aa:.0f}" == f"{max_aa:.0f}":
            ident_string += f"{min_aa:.0f}% aa"
        else:
            ident_string += f"{min_aa:.0f}-{max_aa:.0f}% aa"
        block_strings.append(ident_string)

        # multiple seqeunces for a single chain
        if chain.startswith("H"):
            multiples = sum([len(b.heavies) > 1 for b in self.pairs])
            block_strings.append(f"multiple heavies: {multiples}")
        else:
            multiples = sum([len(b.lights) > 1 for b in self.pairs])
            block_strings.append(f"multiple lights: {multiples}")

        return "\n".join(block_strings)

    def _cdr_alignment_block(self, sequences, dot_alignment=True, color=True):
        aln_data = []

        # get a representative sequence for parsing junction germline info
        notnone_sequences = [s for s in sequences if s is not None]
        ref = notnone_sequences[0]

        # get data for germline line
        d = {"name": "germ", "order": 0}

        # CDR1 germline
        cdr1_aln = semiglobal_alignment(
            ref["cdr1_aa"],
            ref["sequence_alignment_aa"],
        )
        cdr1_aa = ref["germline_alignment_aa"][
            cdr1_aln.target_begin : cdr1_aln.target_end + 1
        ]

        # CDR2 germline
        cdr2_aln = semiglobal_alignment(
            ref["cdr2_aa"],
            ref["sequence_alignment_aa"],
        )
        cdr2_aa = ref["germline_alignment_aa"][
            cdr2_aln.target_begin : cdr2_aln.target_end + 1
        ]

        # v-gene portion of junction
        v_aln = semiglobal_alignment(
            ref["v_sequence"], ref["junction"], gap_open=50, gap_extend=50
        )
        junction_v_length = v_aln.query_end - v_aln.query_begin + 1
        junction_v_nt = ref["v_germline"][-junction_v_length:]
        junction_v = translate(junction_v_nt)

        # j-gene portion of junction
        j_aln = semiglobal_alignment(
            ref["j_sequence"], ref["junction"], gap_open=50, gap_extend=50
        )
        junction_j_length = j_aln.query_end + 1
        junction_j_nt = ref["j_germline"][:junction_j_length]
        junction_j = translate(junction_j_nt, frame=(len(junction_j_nt) % 3) + 1)
        junction_x = "X" * (len(ref["junction_aa"]) - len(junction_v) - len(junction_j))

        # add the germline data
        d["cdr1"] = cdr1_aa
        d["cdr2"] = cdr2_aa
        d["junction"] = junction_v + junction_x + junction_j
        d["junction_v"] = junction_v
        d["junction_j"] = junction_j
        aln_data.append(d)

        # get data for each sequence line
        for i, seq in enumerate(sequences, 1):
            d = {"name": i, "order": i}
            for region in ["cdr1", "cdr2", "junction"]:
                if seq is not None:
                    d[region] = seq[f"{region}_aa"]
                else:
                    d[region] = "X"
            aln_data.append(d)

        # make a dataframe
        df = pd.DataFrame(aln_data).sort_values(by="order")

        # do the alignments
        names = df["name"]
        for region in ["cdr1", "cdr2", "junction"]:
            aln_seqs = [Sequence(s, id=n) for s, n in zip(df[region], names)]
            aln = muscle(aln_seqs)
            # aln_dict = {rec.id: str(rec.seq) for rec in aln}
            aln_dict = {}
            for rec in aln:
                if not str(rec.sequence).replace("-", "").replace("X", ""):
                    aln_dict[rec.id] = " " * len(str(rec.sequence))
                else:
                    aln_dict[rec.id] = str(rec.sequence)
            df[f"{region}_aligned"] = [aln_dict[str(n)] for n in names]

        # make a dot alignment
        for region in ["cdr1", "cdr2"]:
            dot_alns = []
            dot_ref = df[f"{region}_aligned"][0]
            dot_alns.append(dot_ref)
            for name, reg in zip(df["name"], df[f"{region}_aligned"]):
                if name == "germ":
                    continue
                dot_aln = ""
                for r1, r2 in zip(dot_ref, reg):
                    if r1 == r2 == "-":
                        dot_aln += "-"
                    elif r1 == r2:
                        dot_aln += "."
                    else:
                        dot_aln += r2
                dot_alns.append(dot_aln)
            df[f"{region}_dot_aligned"] = dot_alns
        df["junction_dot_aligned"] = df["junction_aligned"]

        # make the block
        block_data = []
        header = []
        sep = " " * self.padding_width
        regions = ["cdr1", "cdr2", "junction"]
        is_dot = "_dot" if dot_alignment else ""
        for _, row in df.iterrows():
            d = []
            if row["name"] == "germ":
                for region in regions:
                    c = self.COLOR_PREFIX[self.REGION_COLORS[region]] if color else ""
                    r = row[f"{region}{is_dot}_aligned"]
                    if region == "junction":
                        vc = (
                            self.COLOR_PREFIX[self.REGION_COLORS["junction_v"]]
                            if color
                            else ""
                        )
                        jc = (
                            self.COLOR_PREFIX[self.REGION_COLORS["junction_j"]]
                            if color
                            else ""
                        )
                        r = r.replace("X", ".")
                        # v = r.split(".")[0]
                        # j = r.split(".")[-1]
                        if "." in r:
                            v = r.split(".")[0]
                            j = r.split(".")[-1]
                        else:
                            v = ""
                            iter_r = iter(r)
                            while len(v.replace("-", "")) < len(row["junction_v"]):
                                v += next(iter_r)
                            j = ""
                            iter_r_rev = iter(r[::-1])
                            while len(j.replace("-", "")) < len(row["junction_j"]):
                                j += next(iter_r_rev)
                            j = j[::-1]
                        n = "." * (len(r) - len(v) - len(j))
                        r = f"{vc}{v}{c}{n}{jc}{j}{c}"
                        # build the junction header
                        vheader = "v"
                        vheader += "-" * (len(v) - 1)
                        # vheader += '>'
                        vheader = vheader[: len(v)]
                        # jheader = '<'
                        jheader = "-" * (len(j) - 1)
                        jheader += "j"
                        jheader = jheader[-len(j) :]
                        header.append(f"{vheader}{' ' * len(n)}{jheader}")
                    else:
                        rname = region if len(region) <= len(r) else region[-1:]
                        h_dashes = max((len(r) - len(rname)) / 2, 0)
                        prefix_dashes = "-" * int(np.floor(h_dashes))
                        suffix_dashes = "-" * int(np.ceil(h_dashes))
                        header.append(f"{prefix_dashes}{rname}{suffix_dashes}")
                        r = f"{c}{r}"
                    d.append(r)
            else:
                d.extend([row[f"{region}{is_dot}_aligned"] for region in regions])
            block_data.append([str(_d) for _d in d])
        # assemble the alignment table
        x = pt.PrettyTable()
        for bd in [header] + block_data:
            x.add_row(bd)
        x.set_style(pt.PLAIN_COLUMNS)
        x.header = False
        x.left_padding_width = self.padding_width
        x.right_padding_width = self.padding_width
        table_string = x.get_string()
        unpadded = "\n".join([l.strip() for l in table_string.split("\n")])
        return unpadded


def fast_tree(alignment, tree, is_aa, show_output=False):
    if is_aa:
        ft_cmd = "fasttree {} > {}".format(alignment, tree)
    else:
        ft_cmd = "fasttree -nt {} > {}".format(alignment, tree)
    ft = sp.Popen(ft_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = ft.communicate()
    if show_output:
        print(ft_cmd)
        print(stdout)
        print(stderr)
    return tree


def group_lineages(pairs, just_pairs=False):
    lineages = {}
    if just_pairs:
        pairs = [p for p in pairs if p.is_pair]
    for p in pairs:
        if p.heavy is not None:
            if "clonify" in p.heavy:
                l = p.heavy["clonify"]["id"]
                if l not in lineages:
                    lineages[l] = []
                lineages[l].append(p)
    return [Lineage(v) for v in lineages.values()]


# class MySequenceItem(QGraphicsRectItem):
#     def __init__(self, seq, seqtype="aa", poswidth=1, posheight=10,
#                  draw_text=False):
#         QGraphicsRectItem.__init__(self)
#         self.seq = seq
#         self.seqtype = seqtype
#         self.poswidth = poswidth
#         self.posheight = posheight
#         if draw_text:
#             self.poswidth = poswidth
#         self.draw_text = draw_text
#         if seqtype == "aa":
#             self.fg = _aafgcolors
#             self.bg = _aabgcolors
#         elif seqtype == "nt":
#             self.fg = _ntfgcolors
#             self.bg = _ntbgcolors
#         self.setRect(0, 0, len(seq) * poswidth, posheight)

#     def paint(self, p, option, widget):
#         x, y = 0, 0
#         qfont = QFont("Courier")
#         current_pixel = 0
#         blackPen = QPen(QColor("black"))
#         for letter in self.seq:
#             if x >= current_pixel:
#                 if self.draw_text and self.poswidth >= 5:
#                     br = QBrush(QColor(self.bg.get(letter, "white")))
#                     p.setPen(blackPen)
#                     p.fillRect(x, 0, self.poswidth, self.posheight, br)
#                     qfont.setPixelSize(min(self.posheight, self.poswidth))
#                     p.setFont(qfont)
#                     p.setBrush(QBrush(QColor("black")))
#                     p.drawText(x, 0, self.poswidth, self.posheight,
#                                Qt.AlignCenter | Qt.AlignVCenter,
#                                letter)
#                 elif letter == "-" or letter == ".":
#                     p.setPen(blackPen)
#                     p.drawLine(x, self.posheight / 2, x + self.poswidth, self.posheight / 2)

#                 else:
#                     br = QBrush(QColor(self.bg.get(letter, "white")))
#                     p.fillRect(x, 0, max(1, self.poswidth), self.posheight, br)
#                     # p.setPen(QPen(QColor(self.bg.get(letter, "black"))))
#                     # p.drawLine(x, 0, x, self.posheight)
#                 current_pixel = int(x)
#             x += self.poswidth

# _aafgcolors = {
#     'A': "#000000",
#     'R': "#000000",
#     'N': "#000000",
#     'D': "#000000",
#     'C': "#000000",
#     'Q': "#000000",
#     'E': "#000000",
#     'G': "#000000",
#     'H': "#000000",
#     'I': "#000000",
#     'L': "#000000",
#     'K': "#000000",
#     'M': "#000000",
#     'F': "#000000",
#     'P': "#000000",
#     'S': "#000000",
#     'T': "#000000",
#     'W': "#000000",
#     'Y': "#000000",
#     'V': "#000000",
#     'B': "#000000",
#     'Z': "#000000",
#     'X': "#000000",
#     '.': "#000000",
#     '-': "#000000",
# }

# _aabgcolors = {
#     'A': "#C8C8C8",
#     'R': "#145AFF",
#     'N': "#00DCDC",
#     'D': "#E60A0A",
#     'C': "#E6E600",
#     'Q': "#00DCDC",
#     'E': "#E60A0A",
#     'G': "#EBEBEB",
#     'H': "#8282D2",
#     'I': "#0F820F",
#     'L': "#0F820F",
#     'K': "#145AFF",
#     'M': "#E6E600",
#     'F': "#3232AA",
#     'P': "#DC9682",
#     'S': "#FA9600",
#     'T': "#FA9600",
#     'W': "#B45AB4",
#     'Y': "#3232AA",
#     'V': "#0F820F",
#     'B': "#FF69B4",
#     'Z': "#FF69B4",
#     'X': "#BEA06E",
#     '.': "#FFFFFF",
#     '-': "#FFFFFF",
# }

_aa_pixel_colors = {
    "A": "#C8C8C8",
    "R": "#145AFF",
    "N": "#00DCDC",
    "D": "#E60A0A",
    "C": "#E6E600",
    "Q": "#00DCDC",
    "E": "#E60A0A",
    "G": "#EBEBEB",
    "H": "#8282D2",
    "I": "#0F820F",
    "L": "#0F820F",
    "K": "#145AFF",
    "M": "#E6E600",
    "F": "#3232AA",
    "P": "#DC9682",
    "S": "#FA9600",
    "T": "#FA9600",
    "W": "#B45AB4",
    "Y": "#3232AA",
    "V": "#0F820F",
    "B": "#FF69B4",
    "Z": "#FF69B4",
    "X": "#BEA06E",
    ".": "#F2F2F2",
    "-": "#FFFFFF",
}

# _ntfgcolors = {
#     'A': '#000000',
#     'G': '#000000',
#     'I': '#000000',
#     'C': '#000000',
#     'T': '#000000',
#     'U': '#000000',
#     '.': "#000000",
#     '-': "#000000",
#     ' ': "#000000"
# }

# _ntbgcolors = {
#     'A': '#A0A0FF',
#     'G': '#FF7070',
#     'I': '#80FFFF',
#     'C': '#FF8C4B',
#     'T': '#A0FFA0',
#     'U': '#FF8080',
#     '.': "#FFFFFF",
#     '-': "#FFFFFF",
#     ' ': "#FFFFFF"
# }

_nt_pixel_colors = {
    "A": "#E12427",
    "C": "#3B7FB6",
    "G": "#63BE7A",
    "T": "#E1E383",
    # 'U': '#E1E383',
    ".": "#F2F2F2",
    "-": "#FFFFFF",
    # ' ': "#FFFFFF"
}
