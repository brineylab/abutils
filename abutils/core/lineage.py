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


from __future__ import absolute_import, division, print_function, unicode_literals

import colorsys
from collections import Counter
import math
import os
import random
import string
import subprocess as sp
import sys

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import ete3

from abstar import run as run_abstar

from .pair import Pair
from .sequence import Sequence
from ..utils.alignment import mafft, muscle
from ..utils.cluster import cluster
from ..utils.color import hex_to_rgb, get_cmap
from ..utils.decorators import lazy_property

# # imports to overload ete2's SequenceItem class
# if sys.version_info[0] > 2:
#     from PyQt5.QtWidgets import QGraphicsRectItem
#     from PyQt5.QtGui import QPen, QColor, QBrush, QFont
#     from PyQt5.QtCore import Qt
# else:
#     from PyQt4.QtGui import QGraphicsRectItem, QPen, QColor, QBrush, QFont
#     from PyQt4.QtCore import Qt


class Lineage(object):
    '''
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

    '''
    def __init__(self, pairs, name=None):
        super(Lineage, self).__init__()
        self.pairs = pairs
        self.rmp_threshold = 0.9
        self.rmp_alt_allowed_mismatches = 1


    def __contains__(self, item):
        if item in self.pair_dict.keys():
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
        '''
        Returns the lineage name, or None if the name cannot be found.
        '''
        clonify_ids = [p.heavy['clonify']['id'] for p in self.heavies if 'clonify' in p.heavy]
        if len(clonify_ids) > 0:
            return clonify_ids[0]
        return None

    @lazy_property
    def just_pairs(self):
        '''
        Returns all lineage Pair objects that contain both
        heavy and light chains.
        '''
        return [p for p in self.pairs if p.is_pair]

    @lazy_property
    def verified_pairs(self):
        '''
        Returns all lineage Pair objects that contain verified pairings.
        '''
        if not hasattr(self.just_pairs[0], 'verified'):
            self.verify_light_chains()
        return [p for p in self.just_pairs if p.verified]

    @lazy_property
    def heavies(self):
        '''
        Returns all lineage Pair objects that contain a heavy chain.
        '''
        return [p for p in self.pairs if p.heavy is not None]

    @lazy_property
    def lights(self):
        '''
        Returns all lineage Pair objects that contain a light chain.
        '''
        return [p for p in self.pairs if p.light is not None]

    @lazy_property
    def uca(self):
        '''
        Calculates the unmutated common ancestor (UCA) of the lineage.

        The UCA is computed by using the germline V(D)J genes as well
        as the N-addition regions of the least mutated sequence.

        If both heavy and light chains exist in the lineage, the UCA
        will be computed for both chains. If not, the UCA will only
        be computed on the chains that exist in the lineage.

        Returns
        -------
        A VaxTools Pair object containing the UCA.
        '''
        return self._calculate_uca()

    @lazy_property
    def rmp(self):
        return self._calculate_recalled_memory_precursor()

    @lazy_property
    def subject(self):
        return self._get_metadata('subject')

    @lazy_property
    def group(self):
        return self._get_metadata('group')

    @lazy_property
    def experiment(self):
        return self._get_metadata('experiment')

    @lazy_property
    def timepoints(self):
        if self.heavies:
            _timepoints = [p.heavy._mongo.get('timepoint', None) for p in self.heavies]
            return list(set([g for g in _timepoints if g is not None]))
        if self.lights:
            _timepoints = [p.light._mongo.get('timepoint', None) for p in self.lights]
            return list(set([g for g in _timepoints if g is not None]))
        return []

    @lazy_property
    def has_insertion(self):
        ins = ['v_ins' in p.heavy for p in self.heavies]
        ins += ['v_ins' in p.light for p in self.lights]
        return True if any(ins) else False

    @lazy_property
    def has_deletion(self):
        dels = ['v_del' in p.heavy for p in self.heavies]
        dels += ['v_del' in p.light for p in self.lights]
        return True if any(dels) else False

    @lazy_property
    def has_indel(self):
        if any([self.has_insertion, self.has_deletion]):
            return True
        return False


    def size(self, pairs_only=False):
        '''
        Calculate the size of the lineage.

        Inputs (optional)
        -----------------
        pairs_only: count only paired sequences

        Returns
        -------
        Lineage size (int)
        '''
        if pairs_only:
            return len(self.just_pairs)
        else:
            return len(self.heavies)


    def verify_light_chains(self, threshold=0.9):
        '''
        Clusters the light chains to identify potentially spurious (non-lineage)
        pairings. Following clustering, all pairs in the largest light chain
        cluster are assumed to be correctly paired. For each of those pairs,
        the <verified> attribute is set to True. For pairs not in the largest
        light chain cluster, the <verified> attribute is set to False.

        Inputs (optional)
        -----------------
        threshold: CD-HIT clustering threshold. Default is 0.9.
        '''
        lseqs = [l.light for l in self.lights]
        clusters = cluster(lseqs, threshold=threshold)
        clusters.sort(key=lambda x: x.size, reverse=True)
        verified_ids = clusters[0].ids
        for p in self.lights:
            p.verified = True if p.name in verified_ids else False


    def dot_alignment(self, seq_field='vdj_nt', name_field='seq_id', uca=None,
            chain='heavy', uca_name='UCA', as_fasta=False, just_alignment=False):
        '''
        Returns a multiple sequence alignment of all lineage sequence with the UCA
        where matches to the UCA are shown as dots and mismatches are shown as the
        mismatched residue.

        Inputs (optional)
        -----------------
        seq_field: the sequence field to be used for alignment. Default is 'vdj_nt'.
        name_field: field used for the sequence name. Default is 'seq_id'.
        chain: either 'heavy' or 'light'. Default is 'heavy'.

        Returns
        -------
        The dot alignment (string)
        '''
        if uca is None:
            uca = self.uca.heavy if chain == 'heavy' else self.uca.light
        uca.id = 'UCA'
        if chain == 'heavy':
            sequences = [p.heavy for p in self.heavies if seq_field in p.heavy]
            if name_field != 'seq_id':
                uca[name_field] = uca['seq_id']
            sequences.append(uca)
            seqs = [(s[name_field], s[seq_field]) for s in sequences]
        else:
            sequences = [p.light for p in self.lights if seq_field in p.light]
            if name_field != 'seq_id':
                uca[name_field] = uca['seq_id']
            sequences.append(uca)
            seqs = [(s[name_field], s[seq_field]) for s in sequences]
        aln = muscle(seqs)
        g_aln = [a for a in aln if a.id == 'UCA'][0]
        dots = [(uca_name, str(g_aln.seq)), ]
        for seq in [a for a in aln if a.id != 'UCA']:
            s_aln = ''
            for g, q in zip(str(g_aln.seq), str(seq.seq)):
                if g == q == '-':
                    s_aln += '-'
                elif g == q:
                    s_aln += '.'
                else:
                    s_aln += q
            dots.append((seq.id, s_aln))
        if just_alignment:
                return [d[1] for d in dots]
        name_len = max([len(d[0]) for d in dots]) + 2
        dot_aln = []
        for d in dots:
            if as_fasta:
                dot_aln.append('>{}\n{}'.format(d[0], d[1]))
            else:
                spaces = name_len - len(d[0])
                dot_aln.append(d[0] + ' ' * spaces + d[1])
        return '\n'.join(dot_aln)


    def pixel(self, seq_field='vdj_nt', figfile=None):
        if 'aa' in seq_field:
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
        plt.imshow(data, cmap=cmap, interpolation='none')
        plt.axis('off')
        if figfile is not None:
            plt.savefig(figfile, bbox_inches='tight', dpi=400)
            plt.close()
        else:
            plt.show()


    # def phylogeny(self, project_dir, aln_file=None, tree_file=None, root=None, seq_field='vdj_nt', aa=False,
    #         root_name=None, show_root_name=False, colors=None, color_function=None, orders=None, order_function=None,
    #         chain='heavy', filter_function=None, just_pairs=False, color_node_labels=False, label_colors=None,
    #         scale=None, branch_vert_margin=None, fontsize=12, show_names=True, name_field='seq_id', show_scale=False,
    #         mirror=False, min_order_fraction=0.1, figname_prefix=None, figname_suffix=None,
    #         linked_alignment=None, alignment_fontsize=11, scale_factor=1, rename_function=None,
    #         alignment_height=50, alignment_width=50, compact_alignment=False, linewidth=1.0,
    #         show_output=False):
    #     '''
    #     Generates a lineage phylogeny figure.

    #     Inputs (required)
    #     -----------------
    #     project_dir: directory for all phylogeny files,
    #         including alignment, tree and figure files

    #     Inputs (optional)
    #     -----------------
    #     aln_file: if a multiple sequence alignment has already been calculated,
    #         passing the path to the alignment file will force Lineage.phylogeny()
    #         to use the supplied msa instead of computing a new one.
    #     tree_file: if a tree file has already been calculated, passing the path
    #         to the pre-computed tree file will force Lineage.phylogeny() to use
    #         the supplied tree file instead of computing a new one.
    #     aa: if True, use amino acid sequences to compute the phylogeny.
    #         Default is False.
    #     root: provide a sequence to be used as the tree root. If not provided,
    #         the UCA will be used to root the tree.
    #     colors: a dictionary with sequence IDs as keys and colors as values. If any
    #         lineage sequences are not in the dict, they will be colored black. If
    #         not provided, all leaves will be colored black. Alternately, if provided
    #         in combination with either <orders> or <order_function>, the dictionary
    #         keys should be orders (integers) instead of sequence IDs.
    #     color_function: provide a function that will be called on each sequence. The
    #         function should accept an AbTools Sequence object and return a color
    #         (as a hex value).
    #     orders: a dictionary with sequence IDs as keys and orders (integers) as values.
    #         If not provided, only the leaf branches will be colored (if <colors> or
    #         <color_function> is provided).
    #     chain: build a phylogeny using the given chain ('heavy' or 'light').
    #         Default is 'heavy'.
    #     filter_function: function used to filter sequences (identity-based clustering, for
    #         example). The function should accept a list of Sequence objects and return
    #         a list of Sequence objects.
    #     just_pairs: if True, compute the phylogeny using only paired sequences.
    #         Default (False) will use all sequences of the appropriate chain, paired or not.
    #     scale: passed to ete3.TreeStyle() to set the scale of the tree figure. Increased
    #         scale results in a wider tree.
    #     branch_vert_margin: passed to ete3.TreeStyle() to set the branch_vertical_margin of
    #         the tree figure. Increased branch_vert_margin results in a taller tree.
    #     fontsize: size of the leaf labels. Default is 12.
    #     show_names: show names of leaf nodes. Options are True (show labels for all leaf nodes),
    #         False (don't show labels for any leaf nodes) or a list of sequence IDs for which
    #         labels should be shown. Default is True.
    #     mirror: flip the orientation of the tree. Default is to draw the tree from left to right.
    #         Setting mirror to True results in the tree being drawn from right to left.
    #     min_order_fraction: minimum fraction of downstream leaves requried to color a branch.
    #         When coloring non-leaf nodes, the earliest 'order' with at least <min_order_fraction>
    #         leaf nodes is used. Default is 0.1 (which corresponds to 10%).
    #     figname_prefix: by default, figures will be named <lineage_id>.pdf. If prefix='prefix_' and
    #         the lineage ID is 'ABC123', the figure file will be named 'prefix_ABC123.pdf'.
    #     figname_suffix: by default, figures will be named <lineage_id>.pdf. If suffix='_suffix' and
    #         the lineage ID is 'ABC123', the figure file will be named 'ABC123_suffix.pdf'.
    #     '''
    #     # seq_field = 'vdj_nt' if seq_field is None else seq_field
    #     project_dir = os.path.abspath(project_dir)
    #     orientation = 1 if mirror else 0
    #     root_name = root_name if root_name is not None else 'root'
    #     # get sequences for phylogeny
    #     if chain == 'heavy':
    #         seq_pool = self.just_pairs if just_pairs else self.heavies
    #         # seqs = [Sequence(p.heavy[seq_field], id=p.heavy[name_field]) for p in seq_pool]
    #         seqs = [p.heavy for p in seq_pool]
    #         if filter_function is not None:
    #             seqs = filter_function(seqs)
    #         if root is None:
    #             root = self.uca.heavy[seq_field]
    #         seqs += [Sequence(root, id=root_name)]
    #     else:
    #         seq_pool = self.just_pairs if just_pairs else self.lights
    #         # seqs = [Sequence(p.light[seq_field], id=p.light[name_field]) for p in seq_pool]
    #         seqs = [p.light for p in seq_pool]
    #         if filter_function is not None:
    #             seqs = filter_function(seqs)
    #         if root is None:
    #             root = self.uca.light[seq_field]
    #         seqs += [Sequence(root, id=root_name)]
    #     # filter sequences
    #     # if filter_function is not None:
    #     #     seqs = filter_function(seqs)
    #     # setup orders
    #     if orders is None:
    #         if order_function is not None:
    #             orders = {seq['seq_id']: order_function(seq) for seq in seqs}
    #     # setup colors
    #     if colors is None:
    #         if color_function is not None:
    #             colors = {seq['seq_id']: color_function(seq) for seq in seqs}
    #         else:
    #             colors = {}
    #     # make msa
    #     if all([aln_file is None, tree_file is None]):
    #         aln_file = os.path.abspath(os.path.join(project_dir, '{}.aln'.format(self.name)))
    #         # muscle(seqs, aln_file, as_file=True)
    #         mafft(seqs, aln_file, as_file=True)
    #     # make treefile
    #     if tree_file is None:
    #         tree_file = os.path.abspath(os.path.join(project_dir, '{}.nw'.format(self.name)))
    #         fast_tree(aln_file, tree_file, is_aa=aa, show_output=show_output)
    #     # make phylogeny
    #     prefix = '' if figname_prefix is None else figname_prefix
    #     suffix = '' if figname_suffix is None else figname_suffix
    #     fig_file = os.path.join(project_dir, '{}{}{}.pdf'.format(prefix, self.name, suffix))
    #     self._make_tree_figure(tree_file, fig_file, colors, orders, root_name, rename_function=rename_function,
    #         show_names=show_names, name_field=name_field, branch_vert_margin=branch_vert_margin, scale=scale,
    #         color_node_labels=color_node_labels, label_colors=label_colors, show_root_name=show_root_name,
    #         tree_orientation=orientation, fontsize=fontsize, min_order_fraction=min_order_fraction,
    #         linked_alignment=linked_alignment, alignment_fontsize=alignment_fontsize, chain=chain,
    #         alignment_height=alignment_height, alignment_width=alignment_width, show_scale=show_scale,
    #         compact_alignment=compact_alignment, scale_factor=scale_factor, linewidth=linewidth)


    def _get_metadata(self, field):
        _metadata = None
        if self.heavies:
            _metadata = [p.heavy._mongo.get(field, None) for p in self.heavies]
            _metadata = [g for g in _metadata if g is not None]
            if len(list(set(_metadata))) == 1:
                metadata = _metadata[0]
            elif len(list(set(_metadata))) > 1:
                mcount = Counter(_metadata)
                metadata = sorted([s for s in mcount.keys()],
                                key=lambda x: mcount[s],
                                reverse=True)[0]
        if self.lights and metadata is None:
            _metadata = [p.light._mongo.get(field, None) for p in self.lights]
            _metadata = [g for g in _metadata if g is not None]
            if len(list(set(_metadata))) == 1:
                metadata = _metadata[0]
            elif len(list(set(_metadata))) > 1:
                mcount = Counter(_metadata)
                metadata = sorted([s for s in mcount.keys()],
                                key=lambda x: mcount[s],
                                reverse=True)[0]
        return metadata


    # def _make_tree_figure(self, tree, fig, colors, orders, root_name, scale=None, branch_vert_margin=None,
    #         fontsize=12, show_names=True, name_field='seq_id', rename_function=None, color_node_labels=False, label_colors=None,
    #         tree_orientation=0, min_order_fraction=0.1, show_root_name=False, chain=None,
    #         linked_alignment=None, alignment_fontsize=11, alignment_height=50, alignment_width=50,
    #         compact_alignment=False, scale_factor=1, linewidth=1, show_scale=False):
    #     if show_names is True:
    #         if chain == 'heavy':
    #             show_names = [p.heavy[name_field] for p in self.pairs if p.heavy is not None]
    #         else:
    #             show_names = [p.light[name_field] for p in self.pairs if p.light is not None]
    #     elif show_names is False:
    #         show_names = []
    #     if show_root_name is True:
    #         show_names.append(root_name)
    #     if linked_alignment is not None:
    #         t = ete3.PhyloTree(tree, alignment=linked_alignment, alg_format='fasta')
    #         ete3.faces.SequenceItem = MySequenceItem
    #     else:
    #         t = ete3.Tree(tree)
    #     t.set_outgroup(t&root_name)
    #     # style the nodes
    #     for node in t.traverse():
    #         if orders is not None:
    #             leaves = node.get_leaf_names()
    #             order_count = Counter([orders[l] for l in leaves])
    #             for order in sorted(order_count.keys()):
    #                 if float(order_count[order]) / len(leaves) >= min_order_fraction:
    #                     color = colors[order]
    #                     break
    #         else:
    #             color = colors.get(node.name, '#000000')
    #         if linked_alignment is not None:
    #             node.add_feature('aln_fontsize', alignment_fontsize)
    #             node.add_feature('aln_height', alignment_height)
    #             node.add_feature('aln_width', alignment_width)
    #             node.add_feature('fontsize', fontsize)
    #             node.add_feature('format', 'seq')
    #             node.add_feature('scale_factor', scale_factor)
    #         style = ete3.NodeStyle()
    #         style['size'] = 0
    #         style['vt_line_width'] = float(linewidth)
    #         style['hz_line_width'] = float(linewidth)
    #         style['vt_line_color'] = color
    #         style['hz_line_color'] = color
    #         style['vt_line_type'] = 0
    #         style['hz_line_type'] = 0
    #         # else:
    #         #     style['size'] = 0
    #         #     style['vt_line_width'] = float(linewidth)
    #         #     style['hz_line_width'] = float(linewidth)
    #         #     style['vt_line_color'] = color
    #         #     style['hz_line_color'] = color
    #         #     style['vt_line_type'] = 0
    #         #     style['hz_line_type'] = 0
    #         if node.name in show_names:
    #             if color_node_labels:
    #                 if label_colors is None:
    #                     node_color = color
    #                 elif type(label_colors) == dict:
    #                     node_color = label_colors.get(node.name, '#000000')
    #                 elif type(label_colors) in [list, tuple]:
    #                     node_color = color if node.name in label_colors else '#000000'
    #                 else:
    #                     node_color = '#000000'
    #             else:
    #                 node_color = '#000000'
    #             node_name = node.name if rename_function is None else rename_function(node.name)
    #             tf = ete3.TextFace(node_name,
    #                                fsize=fontsize,
    #                                fgcolor=node_color)
    #             # tf.fsize = fontsize
    #             node.add_face(tf, column=0)
    #             # style['fgcolor'] = hex_to_rgb(node_color)
    #         # else:
    #         #     if hasattr(node, "sequence"):
    #         #         node.add_face(ete3.SeqMotifFace(seq=node.sequence,
    #         #                                         seqtype="aa",
    #         #                                         height=50,
    #         #                                         seq_format="seq"), column=0, position="aligned")
    #         node.set_style(style)
    #     t.dist = 0
    #     ts = ete3.TreeStyle()
    #     if linked_alignment is not None:
    #         ts.layout_fn = self._phyloalignment_layout_function
    #     ts.orientation = tree_orientation
    #     ts.show_leaf_name = False
    #     if scale is not None:
    #         ts.scale = int(scale)
    #     if branch_vert_margin is not None:
    #         ts.branch_vertical_margin = float(branch_vert_margin)
    #     ts.show_scale = show_scale
    #     # ladderize
    #     t.ladderize()
    #     # render the tree
    #     t.render(fig, tree_style=ts)


    # def _phyloalignment_layout_function(self, node):
    #     leaf_color = "#000000"
    #     node.img_style["shape"] = "circle"
    #     if hasattr(node, "evoltype"):
    #         if node.evoltype == 'D':
    #             node.img_style["fgcolor"] = "#FF0000"
    #             node.img_style["hz_line_color"] = "#FF0000"
    #             node.img_style["vt_line_color"] = "#FF0000"
    #         elif node.evoltype == 'S':
    #             node.img_style["fgcolor"] = "#1d176e"
    #             node.img_style["hz_line_color"] = "#1d176e"
    #             node.img_style["vt_line_color"] = "#1d176e"
    #         elif node.evoltype == 'L':
    #             node.img_style["fgcolor"] = "#777777"
    #             node.img_style["vt_line_color"] = "#777777"
    #             node.img_style["hz_line_color"] = "#777777"
    #             node.img_style["hz_line_type"] = 1
    #             node.img_style["vt_line_type"] = 1
    #             leaf_color = "#777777"

    #     if node.is_leaf():
    #         node.img_style["shape"] = "square"
    #         node.img_style["size"] = 0
    #         if hasattr(node, "sequence"):
    #             if node.name == 'root':
    #                 bg_colors, fg_colors = self._get_phyloalignment_colors(root=True)
    #                 node.img_style["fgcolor"] = '#d3d3d3'
    #                 SequenceFace = ete3.faces.SeqMotifFace(node.sequence, seqtype="aa", seq_format='seq',
    #                     height=node.aln_height, width=node.aln_width, scale_factor=node.scale_factor)
    #                 ete3.faces.add_face_to_node(SequenceFace, node, 1, aligned=True)
    #                 node.name = ' UCA  '
    #                 ete3.faces.add_face_to_node(ete3.faces.AttrFace("name", "Arial", node.fontsize, '#000000', None),
    #                                             node, 0)
    #             else:
    #                 bg_colors, fg_colors = self._get_phyloalignment_colors()
    #                 node.img_style["fgcolor"] = '#000000'
    #                 SequenceFace = ete3.faces.SeqMotifFace(node.sequence, seqtype="aa", seq_format='seq',
    #                     height=node.aln_height, width=node.aln_width, scale_factor=node.scale_factor)
    #                 ete3.faces.add_face_to_node(SequenceFace, node, 1, aligned=True)
    #     else:
    #         node.img_style["size"] = 0


    # def _get_phyloalignment_colors(self, root=False):
    #     bg = '#000000'
    #     fg = '#FFFFFF'
    #     bg_colors = {c: bg for c in string.ascii_uppercase}
    #     bg_colors['.'] = '#FFFFFF'
    #     bg_colors['-'] = '#d3d3d3'
    #     fg_colors = {c: fg for c in string.ascii_uppercase}
    #     fg_colors['.'] = '#000000'
    #     fg_colors['-'] = '#000000'
    #     return bg_colors, fg_colors


    def _calculate_uca(self, paired_only=False):
        if paired_only:
            heavies = self.just_pairs
            lights = self.just_pairs
        else:
            heavies = self.heavies
            lights = self.lights
        # heavy chain UCA
        if len(heavies) >= 1:
            lmhc = sorted(heavies, key=lambda x: x.heavy['nt_identity']['v'], reverse=True)[0].heavy
            hc_uca = run_abstar(('UCA', lmhc['vdj_germ_nt']), isotype=False)
        else:
            hc_uca = None
        # light chain UCA
        if len(lights) >= 1:
            lmlc = sorted(lights, key=lambda x: x.light['nt_identity']['v'], reverse=True)[0].light
            lc_uca = run_abstar(('UCA', lmlc['vdj_germ_nt']), isotype=False)
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
        rmp = ''
        seqs = [(s['seq_id'], s['vdj_nt']) for s in sequences]
        aln = muscle(seqs)
        g_aln = [a for a in aln if a.id == 'UCA'][0]
        query_seqs = [str(a.seq) for a in aln if a.id != 'UCA']
        for i, g in enumerate(g_aln):
            qcounts = Counter([q[i] for q in query_seqs])
            qmax = sorted(qcounts.keys(),
                          key=lambda x: qcounts[x],
                          reverse=True)[0]
            qmax_fraction = float(qcounts[qmax]) / sum(qcounts.values())
            qmax_alt_mismatches = sum(qcounts.values()) - qcounts[qmax]
            if any([qmax_fraction >= self.rmp_threshold,
                    qmax_alt_mismatches <= self.rmp_alt_allowed_mismatches]):
                rmp += qmax
            else:
                rmp += g
        return run_abstar(Sequence(rmp.replace('-', ''), id='RMP'))


    def _germline_field_map(self, field):
        fmap = {'fr1_nt': 'fr1_germ_nt',
                'fr2_nt': 'fr2_germ_nt',
                'fr3_nt': 'fr3_germ_nt',
                'fr4_nt': 'fr4_germ_nt',
                'cdr1_nt': 'cdr1_germ_nt',
                'cdr2_nt': 'cdr2_germ_nt',
                'fr1_aa': 'fr1_germ_aa',
                'fr2_aa': 'fr2_germ_aa',
                'fr3_aa': 'fr3_germ_aa',
                'fr4_aa': 'fr4_germ_aa',
                'cdr1_aa': 'cdr1_germ_aa',
                'cdr2_aa': 'cdr2_germ_aa',
                'vdj_nt': 'vdj_germ_nt',
                'vdj_aa': 'vdj_germ_aa'}
        return fmap.get(field.lower(), None)


def fast_tree(alignment, tree, is_aa, show_output=False):
    if is_aa:
        ft_cmd = 'fasttree {} > {}'.format(alignment, tree)
    else:
        ft_cmd = 'fasttree -nt {} > {}'.format(alignment, tree)
    ft = sp.Popen(ft_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = ft.communicate()
    if show_output:
        print(ft_cmd)
        print(stdout)
        print(stderr)
    return tree


def donut(lineages, figfile=None, figsize=(6, 6), pairs_only=False, monochrome_color=None, singleton_color='lightgrey', shuffle_colors=False, seed=1234,
          text_kws={}, pie_kws={}, fontsize=28):
    lineages = sorted(lineages, key=lambda x: x.size(pairs_only), reverse=True)
    non_singletons = [l for l in lineages if l.size(pairs_only) > 1]
    singleton_count = sum([1 for l in lineages if l.size(pairs_only) == 1])
    lineage_sizes = [l.size(pairs_only) for l in lineages if l.size(pairs_only) > 1] + [singleton_count]
    if monochrome_color is not None:
        colors = _get_monochrome_colors(monochrome_color, len(non_singletons))
        # we shuffle the colors differently if we're using a monochrome palette, because
        # we want to have the first color (largest lineage) always be the user-supplied
        # monochrome_color. We only want to shuffle the colors starting with the second one.
        if shuffle_colors:
            primary = colors[0]
            secondary = colors[1:]
            random.seed(seed)
            random.shuffle(secondary)
            colors = [primary] + secondary
    else:
        colors = _get_donut_colors(len(non_singletons))
        if shuffle_colors:
            random.seed(seed)
            random.shuffle(colors)
    colors += [singleton_color]

    # _colors = _get_donut_colors(lineages, len(lineage_sizes), color, singleton_color, shuffle_colors, seed)
    # for l, c in zip(lineages, _colors):
    #     l.color = c
    # colors = [l.color for l in lineages if l.size(pairs_only) > 1] + [singleton_color]

    fig = plt.figure(figsize=figsize)
    ax = plt.gca()
    # fig, ax = plt.subplots()
    ax.axis('equal')
    width = 0.55
    kwargs = dict(colors=colors, startangle=90)
    for k, v in pie_kws.items():
        kwargs[k] = v
    inside, _ = ax.pie(lineage_sizes, radius=1, pctdistance=1 - width / 2, **kwargs)
    plt.setp(inside, width=width, edgecolor='white')

    for w in inside:
        w.set_linewidth(2)

    kwargs = dict(size=fontsize, color='k', va='center', fontweight='bold')
    for k, v in text_kws.items():
        kwargs[k] = v
    ax.text(0, 0, str(sum(lineage_sizes)), ha='center', **kwargs)

    plt.tight_layout()

    if figfile is not None:
        plt.savefig(figfile)
    else:
        plt.show()


def _get_donut_colors(N):
    HSV_tuples = [(x * 1.0 / N, 0.8, 0.9) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)[::-1]
    return RGB_tuples


def _get_monochrome_colors(monochrome_color, N):
    cmap = get_cmap(from_color=monochrome_color)
    # this is a bit convoluted, but what's happening is we're getting different colormap
    # values (which range from 0 to 1). Calling cmap(i) returns an rgba tuple, but we just need
    # the rbg, so we drop the a. To make sure that one of the colors isn't pure white,
    # we ask np.linspace() for one more value than we need, reverse the list of RGB tuples
    # so that it goes from dark to light, and drop the lightest value
    RGB_tuples = [cmap(i)[:-1] for i in np.linspace(0, 1, N + 1)][::-1][:-1]
    return RGB_tuples



def group_lineages(pairs, just_pairs=False):
    lineages = {}
    if just_pairs:
        pairs = [p for p in pairs if p.is_pair]
    for p in pairs:
        if p.heavy is not None:
            if 'clonify' in p.heavy:
                l = p.heavy['clonify']['id']
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

# _aa_pixel_colors = {
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
#     '.': "#F2F2F2",
#     '-': "#FFFFFF",
# }

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

# _nt_pixel_colors = {
#     'A': '#E12427',
#     'C': '#3B7FB6',
#     'G': '#63BE7A',
#     'T': '#E1E383',
#     # 'U': '#E1E383',
#     '.': "#F2F2F2",
#     '-': "#FFFFFF",
#     # ' ': "#FFFFFF"
# }
