#!/usr/bin/env python
# filename: phylogeny.py


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
from copy import copy, deepcopy
import math
import os
import random
import shutil
import string
import subprocess as sp
import sys

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import ete3

from Bio import AlignIO, Phylo

import abstar
from abstar.core.germline import get_imgt_germlines

from ..core.pair import Pair
from ..core.sequence import Sequence
from .alignment import mafft, muscle
from .cluster import cluster
from .color import hex_to_rgb, get_cmap
from .decorators import lazy_property

# # imports to overload ete3's SequenceItem class
# if sys.version_info[0] > 2:
#     from PyQt5.QtWidgets import QGraphicsRectItem
#     from PyQt5.QtGui import QPen, QColor, QBrush, QFont
#     from PyQt5.QtCore import Qt
# else:
#     from PyQt4.QtGui import QGraphicsRectItem, QPen, QColor, QBrush, QFont
#     from PyQt4.QtCore import Qt


if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
else:
    STR_TYPES = [str, unicode]


def phylogeny(sequences=None, project_dir=None, name=None, aln_file=None, tree_file=None,
        seq_field=None, name_field=None, aa=False, species='human', unrooted=False, ladderize=True,
        root=None, root_name=None, show_root_name=False, color_dict=None, color_function=None,
        order_dict=None, order_function=None, color_node_labels=False, label_colors=None,
        scale=None, branch_vert_margin=None, fontsize=12, show_names=True, show_scale=False,
        mirror=False, min_order_fraction=0.1, figname_prefix=None, figname_suffix=None,
        # linked_alignment=None, alignment_fontsize=11, alignment_height=50, alignment_width=50, 
        compact_alignment=False, scale_factor=1, rename_function=None, linewidth=1.0,
        delete_nodes=None, quiet=True):
    '''
    Generates a lineage phylogeny figure.

    Args:

        sequences (list(Sequence)): A list of ``Sequence`` objects from which a phylogeny
            will be calculated. Strictly speaking, they do not need to be ``Sequence`` objects,
            rather, any object that contains the sequence name as the ``id`` attribute (or
            by dictionary-style lookup using the provided ``name_field``) and contains the
            sequence as the ``sequence`` attribute (or by dictionary-stype lookup using the
            provided ``seq_field``).

        project_dir (str): directory into which all phylogeny files will be deposited,
            including alignment, tree and figure files.

        name (str): Name to be used for naming alignment, tree, and phylogeny files. If not
            provided, a random name will be generated.

        aln_file (str): if a multiple sequence alignment has already been calculated,
            passing the path to the alignment file (in FASTA format) will force Lineage.phylogeny()
            to use the supplied msa instead of computing a new one.

        tree_file (str): if a tree file has already been calculated, passing the path
            to the pre-computed tree file will force ``phylogeny()`` to use
            the supplied tree file instead of computing a new one. It is important to note that
            only sequence names will be parsed from the tree_file, so if ``order_function`` or
            ``color_function`` is also provided, ensure that these functions only require the
            sequence ID rather than the entire sequence.

        aa (bool): if True, use amino acid sequences to compute the phylogeny.
            Default is False.

        root (Sequence, str: The root can be provided either as a ``Sequence`` object (if ``sequences``
            are being provided) or as the name of a sequence that can be found either in
            ``sequences`` or in the provided ``aln_file`` or ``tree_file``. Note that if
            either ``aln_file`` or ``tree_file`` are provided, the root must be provided
            as the sequence name, not as a ``Sequence`` object (as the root sequence must
            already be included in either ``aln_file`` or ``tree_file``. If the root is not
            provided, the germline V-gene sequence of the

        color_dict (dict): Dictionary with sequence IDs as keys and colors (hex format) as values. If any
            sequence IDs are not found in the dict, they will be colored black. If neither ``color_dict`` nor
            ``color_function`` is provided, all leaves will be colored black.

        color_function (func): Function that that accepts a ``Sequence`` object and returns the color
            (as a hex code). If ``color_dict`` is also provided, ``color_function`` is ignored. Additionally,
            ``color_function`` will only be used if ``sequences`` are provided. If ``sequences`` are not provided
            (instead using ``aln_file` or ``tree_file``), ``color_dict`` must be used instead of ``color_function``.

        orders: a dictionary with sequence IDs as keys and orders (integers) as values.
            If not provided, only the leaf branches will be colored (if <colors> or
            <color_function> is provided).

        chain: build a phylogeny using the given chain ('heavy' or 'light').
            Default is 'heavy'.

        filter_function: function used to filter sequences (identity-based clustering, for
            example). The function should accept a list of Sequence objects and return
            a list of Sequence objects.

        just_pairs: if True, compute the phylogeny using only paired sequences.
            Default (False) will use all sequences of the appropriate chain, paired or not.

        scale (float): passed to ete3.TreeStyle() to set the scale of the tree figure. Increased
            scale results in a wider tree.

        branch_vert_margin (int): passed to ete3.TreeStyle() to set the branch_vertical_margin of
            the tree figure. Increased branch_vert_margin results in a taller tree.

        fontsize: size of the leaf labels. Default is 12.

        show_names: show names of leaf nodes. Options are True (show labels for all leaf nodes),
            False (don't show labels for any leaf nodes) or a list of sequence IDs for which
            labels should be shown. Default is True.

        mirror: flip the orientation of the tree. Default is to draw the tree from left to right.
            Setting mirror to True results in the tree being drawn from right to left.

        min_order_fraction: minimum fraction of downstream leaves requried to color a branch.
            When coloring non-leaf nodes, the earliest 'order' with at least <min_order_fraction>
            leaf nodes is used. Default is 0.1 (which corresponds to 10%).

        figname_prefix: by default, figures will be named <lineage_id>.pdf. If prefix='prefix_' and
            the lineage ID is 'ABC123', the figure file will be named 'prefix_ABC123.pdf'.

        figname_suffix: by default, figures will be named <lineage_id>.pdf. If suffix='_suffix' and
            the lineage ID is 'ABC123', the figure file will be named 'ABC123_suffix.pdf'.
    '''

    if project_dir is None:
        print('\nERROR: project_dir is required\n')
        sys.exit(1)
    else:
        project_dir = os.path.abspath(project_dir)

    # make a name if one isn't provided
    if name is None:
        name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

    # if sequences are provided, need to process them
    if sequences is not None and all([arg is None for arg in [aln_file, tree_file]]):
        sequences = deepcopy(sequences)
        root = copy(root)

        # if custom seq_field is specified, copy to the .sequence attribute
        if seq_field is not None:
            if not all([seq_field in list(s.annotations.keys()) for s in sequences]):
                print('\nERROR: {} is not present in all of the supplied sequences.\n'.format(seq_field))
                sys.exit(1)
            for s in sequences:
                s.alignment_sequence = s[seq_field]
        else:
            for s in sequences:
                s.alignment_sequence = s.sequence

        # if custom name_field is specified, copy to the .id attribute
        if name_field is not None:
            if not all([name_field in list(s.annotations.keys()) for s in sequences]):
                print('\nERROR: {} is not present in all of the supplied sequences.\n'.format(name_field))
                sys.exit(1)
            for s in sequences:
                s.alignment_id = s[name_field]
        else:
            for s in sequences:
                s.alignment_id = s.id

        # parse the root sequence
        if unrooted:
            root = None
            root_name = None
        elif root is None:
            if not quiet:
                print('\nRoot sequence was was not provided. Using the germline V-gene.')
            if not all(['v_gene' in list(s.annotations.keys()) for s in sequences]):
                print('\nInput sequences to not appear to be AbStar annotated. Annotating now...')
                sequences = abstar.run(*[(s.id, s.sequence) for s in sequences])
                print('Done.')
            if not all(['full' in list(s['v_gene'].keys()) for s in sequences]):
                print('\nInput sequences to not appear to be AbStar annotated. Annotating now...')
                sequences = abstar.run(*[(s.id, s.sequence) for s in sequences])
                print('Done.')
            top_vgene = sorted(list(Counter([s['v_gene']['full'] for s in sequences]).items()),
                               key=lambda x: x[1],
                               reverse=True)[0][0]
            vgene = get_imgt_germlines(species, 'V', gene=top_vgene)
            if aa:
                root = Sequence(vgene.ungapped_aa_sequence, id=top_vgene)
            else:
                root = Sequence(vgene.ungapped_nt_sequence, id=top_vgene)
            root.alignment_id = root.id
            root.alignment_sequence = root.sequence
            if not quiet:
                print('Top V-gene: {}'.format(root.alignment_id))
                print(root.alignment_sequence)
        elif type(root) in STR_TYPES:
            root = [s for s in sequences if s.alignment_id == root][0]
            if not root:
                print('\nERROR: The name of the root sequence ({}) was not found in the list of input sequences.'.format(root))
                print('\n')
                sys.exit(1)
            sequences = [s for s in sequences if s.alignment_id != root.alignment_id]
        elif type(root) == Sequence:
            if seq_field is not None:
                if seq_field not in list(root.anotations.keys()):
                    print('\nERROR: {} is not present in the supplied root sequence.\n'.format(seq_field))
                    sys.exit(1)
                root.alignment_sequence = root[seq_field]
            if name_field is not None:
                if name_field not in list(root.anotations.keys()):
                    print('\nERROR: {} is not present in the supplied root sequence.\n'.format(name_field))
                    sys.exit(1)
                root.alignment_id = root[name_field]
            sequences = [s for s in sequences if s.alignment_id != root.alignment_id]
        else:
            print('\nERROR: If root is provided, it must be the name of a sequence \
            found in the supplied list of sequences or it must be a Sequence object.')
            print('\n')
            sys.exit(1)
        if not unrooted:
            if root_name is not None:
                root.alignment_id = root_name
            else:
                root_name = root.alignment_id
            sequences.append(root)

    # parse sequences from aln_file, if provided
    elif aln_file is not None:
        if not unrooted and type(root) not in STR_TYPES:
            print('\nERROR: If providing an aln_file, the name of the root sequence must \
            be provided (as a string) using the root keyword argument')
            print('\n')
            sys.exit(1)
        _sequences = []
        _root = None
        for rec in AlignIO.read(open(aln_file), 'fasta'):
            s = str(rec.seq).replace('-', '')
            if rec.id == root:
                _root = Sequence(s, rec.id)
                _root.alignment_id = _root.id
            else:
                _s = Sequence(s, id=rec.id)
                _s.alignment_id = rec.id
                _sequences.append(_s)
        if sequences is None:
            sequences = _sequences
        else:
            sequence_ids = [s.id for s in sequences]
            if any([_s.alignment_id not in sequence_ids for _s in _sequences]):
                print('\nWARNING: Sequences were found in the alignment file that were not included \
                      in the input sequence list. This may cause problems.')
            for s in sequences:
                s.alignment_id = s.id
                s.alignment_sequence = s.sequence
        if unrooted:
            root = None
            root_name = None
        else:
            if _root is None:
                print('\nERROR: The specified root ({}) was not found in the provided alignment file.'.format(root))
                print('\n')
                sys.exit(1)
            root = _root
            if root_name is not None:
                root.alignment_id = root_name
            else:
                root_name = root.alignment_id
            sequences = [s for s in sequences if all([s.alignment_id != name for name in [root.id, root.alignment_id]])]
            sequences.append(root)

    # parse sequences from tree_file, if provided
    elif tree_file is not None:
        if not unrooted and type(root) not in STR_TYPES:
            print('\nERROR: If providing a tree_file, the name of the root sequence must \
            be provided (as a string) using the root keyword argument')
            print('\n')
            sys.exit(1)
        _sequences = []
        _root = None
        tree = Phylo.read(open(tree_file), 'newick')
        for leaf in tree.get_terminals():
            s = ''
            if leaf.name == root:
                _root = Sequence(s, leaf.name)
                _root.alignment_id = _root.id
            else:
                _s = Sequence(s, id=leaf.name)
                _s.alignment_id = leaf.name
                _sequences.append(_s)
        if sequences is None:
            sequences = _sequences
        else:
            sequence_ids = [s.id for s in sequences]
            if any([_s.alignment_id not in sequence_ids for _s in _sequences]):
                print('\nWARNING: Sequences were found in the alignment file that were not included \
                      in the input sequence list. This may cause problems.')
            for s in sequences:
                s.alignment_id = s.id
                s.alignment_sequence = s.sequence
        if unrooted:
            root = None
            root_name = None
        elif _root is None:
            print('\nERROR: The specified root ({}) was not found in the provided tree file.'.format(root))
            print('\n')
            sys.exit(1)
        else:
            root = _root
            if root_name is not None:
                root.alignment_id = root_name
            else:
                root_name = root.alignment_id
            sequences = [s for s in sequences if all([s.alignment_id != name for name in [root.id, root.alignment_id]])]
            sequences.append(root)

    # set up colors and color ordering
    if order_dict is None:
        if order_function is not None:
            order_dict = {seq.alignment_id: order_function(seq) for seq in sequences}
    if color_dict is None:
        if color_function is not None:
            color_dict = {seq.alignment_id: color_function(seq) for seq in sequences}
    if color_dict is None:
        color_dict = {}

    # make msa (if necessary)
    if all([aln_file is None, tree_file is None]):
        aln_file = os.path.abspath(os.path.join(project_dir, '{}.aln'.format(name)))
        # muscle(seqs, aln_file, as_file=True)
        do_print = False if quiet else True
        if do_print:
            print('\n')
        seqs = [(s.alignment_id, s.alignment_sequence) for s in sequences]
        mafft(seqs, aln_file, as_file=True, print_stdout=do_print, print_stderr=do_print)

    # make treefile (if necessary)
    if tree_file is None:
        tree_file = os.path.abspath(os.path.join(project_dir, '{}.nw'.format(name)))
        fasttree(aln_file, tree_file, is_aa=aa, quiet=quiet)

    # make phylogeny
    prefix = '' if figname_prefix is None else figname_prefix
    suffix = '' if figname_suffix is None else figname_suffix
    fig_file = os.path.join(project_dir, '{}{}{}.pdf'.format(prefix, name, suffix))
    _make_tree_figure(tree_file,
                      fig_file,
                      color_dict,
                      order_dict,
                      None if root is None else root.alignment_id,
                      rename_function=rename_function,
                      show_names=show_names,
                      name_field=name_field,
                      branch_vert_margin=branch_vert_margin,
                      scale=scale,
                      color_node_labels=color_node_labels,
                      label_colors=label_colors,
                      show_root_name=show_root_name,
                      tree_orientation=1 if mirror else 0,
                      fontsize=fontsize,
                      min_order_fraction=min_order_fraction,
                    #   linked_alignment=linked_alignment,
                    #   alignment_fontsize=alignment_fontsize,
                    #   alignment_height=alignment_height,
                    #   alignment_width=alignment_width,
                      show_scale=show_scale,
                      compact_alignment=compact_alignment,
                      scale_factor=scale_factor,
                      linewidth=linewidth,
                      ladderize=ladderize,
                      delete_nodes=delete_nodes)


def fasttree(alignment, tree_file, is_aa=False, quiet=True):
    if is_aa:
        ft_cmd = 'fasttree {} > {}'.format(alignment, tree_file)
    else:
        ft_cmd = 'fasttree -nt {} > {}'.format(alignment, tree_file)
    ft = sp.Popen(ft_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = ft.communicate()
    if not quiet:
        print(ft_cmd)
        print(stdout)
        print(stderr)
    return tree_file


def lsd(tree, output_file=None, dates_file=None, outgroup_file=None,
        with_constraints=True, with_weights=True, reestimate_root_position=None, quiet=True):
    lsd_cmd = 'lsd -i {}'.format(os.path.abspath(tree))
    if output_file is not None:
        lsd_cmd += ' -o {}'.format(os.path.abspath(output_file))
    if dates_file is not None:
        lsd_cmd += ' -d {}'.format(os.path.abspath(dates_file))
    if outgroup_file is not None:
        lsd_cmd += ' -g {}'.format(os.path.abspath(outgroup_file))
    if with_constraints:
        lsd_cmd += ' -c'
    if with_weights:
        lsd_cmd += ' -v'
    if reestimate_root_position is not None:
        lsd_cmd += ' -r {}'.format(reestimate_root_position)
    p = sp.Popen(lsd_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if not quiet:
        print(lsd_cmd)
        print(stdout)
        print(stderr)
    return output_file


def igphyml(input_file=None, tree_file=None, root=None, verbose=False):
    '''
    Computes a phylogenetic tree using IgPhyML.

    .. note::
        
        IgPhyML must be installed. It can be downloaded from https://github.com/kbhoehn/IgPhyML.

    Args:

        input_file (str): Path to a Phylip-formatted multiple sequence alignment. Required.

        tree_file (str): Path to the output tree file.

        root (str): Name of the root sequence. Required.

        verbose (bool): If `True`, prints the standard output and standard error for each IgPhyML run. 
            Default is `False`.
    '''

    if shutil.which('igphyml') is None:
        raise RuntimeError('It appears that IgPhyML is not installed.\nPlease install and try again.')
    
    # first, tree topology is estimated with the M0/GY94 model
    igphyml_cmd1 = 'igphyml -i {} -m GY -w M0 -t e --run_id gy94'.format(aln_file)
    p1 = sp.Popen(igphyml_cmd1, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout1, stderr1 = p1.communicate()
    if verbose:
        print(stdout1 + '\n')
        print(stderr1 + '\n\n')
    intermediate = input_file + '_igphyml_tree.txt_gy94'

    # now  we fit the HLP17 model once the tree topology is fixed
    igphyml_cmd2 = 'igphyml -i {0} -m HLP17 --root {1} -o lr -u {}_igphyml_tree.txt_gy94 -o {}'.format(input_file,
                                                                                                       root,
                                                                                                       tree_file)
    p2 = sp.Popen(igphyml_cmd2, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout2, stderr2 = p2.communicate()
    if verbose:
        print(stdout2 + '\n')
        print(stderr2 + '\n')
    return tree_file + '_igphyml_tree.txt'


def _make_tree_figure(tree, fig, colors, orders, root_name, scale=None, branch_vert_margin=None,
        fontsize=12, show_names=True, name_field='seq_id', rename_function=None, color_node_labels=False, label_colors=None,
        tree_orientation=0, min_order_fraction=0.1, show_root_name=False, chain=None,
        # linked_alignment=None, alignment_fontsize=11, alignment_height=50, alignment_width=50,
        compact_alignment=False, scale_factor=1, linewidth=1, show_scale=False, ladderize=True, delete_nodes=None):
    if delete_nodes is None:
        delete_nodes = []
    elif type(delete_nodes) in STR_TYPES:
        delete_nodes = [delete_nodes, ]
    if show_root_name is True:
        show_names.append(root_name)
    # if linked_alignment is not None:
    #     t = ete3.PhyloTree(tree, alignment=linked_alignment, alg_format='fasta')
    #     ete3.faces.SequenceItem = MySequenceItem
    else:
        t = ete3.Tree(tree)
    if root_name is not None:
        t.set_outgroup(t&root_name)
    # style the nodes
    for node in t.traverse():
        if node.name in delete_nodes:
            node.delete()
            continue
        if orders is not None:
            leaves = node.get_leaf_names()
            order_count = Counter([orders[l] for l in leaves])
            for order in sorted(order_count.keys()):
                if float(order_count[order]) / len(leaves) >= min_order_fraction:
                    color = colors[order]
                    break
        else:
            color = colors.get(node.name, '#000000')
        if linked_alignment is not None:
            node.add_feature('aln_fontsize', alignment_fontsize)
            node.add_feature('aln_height', alignment_height)
            node.add_feature('aln_width', alignment_width)
            node.add_feature('fontsize', fontsize)
            node.add_feature('format', 'seq')
            node.add_feature('scale_factor', scale_factor)
        style = ete3.NodeStyle()
        style['size'] = 0
        style['vt_line_width'] = float(linewidth)
        style['hz_line_width'] = float(linewidth)
        style['vt_line_color'] = color
        style['hz_line_color'] = color
        style['vt_line_type'] = 0
        style['hz_line_type'] = 0
        if show_names is True:
            tf = _build_node_text_face(node, color_node_labels, color, label_colors, fontsize, rename_function)
            node.add_face(tf, column=0)
        elif node.name in show_names:
            tf = _build_node_text_face(node, color_node_labels, color, label_colors, fontsize, rename_function)
            node.add_face(tf, column=0)
        node.set_style(style)
    t.dist = 0
    ts = ete3.TreeStyle()
    # if linked_alignment is not None:
    #     ts.layout_fn = _phyloalignment_layout_function
    ts.orientation = tree_orientation
    ts.show_leaf_name = False
    if scale is not None:
        ts.scale = int(scale)
    if branch_vert_margin is not None:
        ts.branch_vertical_margin = float(branch_vert_margin)
    ts.show_scale = show_scale
    if ladderize:
        t.ladderize()
    t.render(fig, tree_style=ts)


def _build_node_text_face(node, color_node_labels, color, label_colors, fontsize, rename_function):
    if color_node_labels:
        if label_colors is None:
            node_color = color
        elif type(label_colors) == dict:
            node_color = label_colors.get(node.name, '#000000')
        elif type(label_colors) in [list, tuple]:
            node_color = color if node.name in label_colors else '#000000'
        else:
            node_color = '#000000'
    else:
        node_color = '#000000'
    node_name = node.name if rename_function is None else rename_function(node.name)
    tf = ete3.TextFace(node_name,
                       fsize=fontsize,
                       fgcolor=node_color)
    return tf


# def _phyloalignment_layout_function(node):
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
#                 bg_colors, fg_colors = _get_phyloalignment_colors(root=True)
#                 node.img_style["fgcolor"] = '#d3d3d3'
#                 SequenceFace = ete3.faces.SeqMotifFace(node.sequence, seqtype="aa", seq_format='seq',
#                     height=node.aln_height, width=node.aln_width, scale_factor=node.scale_factor)
#                 ete3.faces.add_face_to_node(SequenceFace, node, 1, aligned=True)
#                 node.name = ' UCA  '
#                 ete3.faces.add_face_to_node(ete3.faces.AttrFace("name", "Arial", node.fontsize, '#000000', None),
#                                             node, 0)
#             else:
#                 bg_colors, fg_colors = _get_phyloalignment_colors()
#                 node.img_style["fgcolor"] = '#000000'
#                 SequenceFace = ete3.faces.SeqMotifFace(node.sequence, seqtype="aa", seq_format='seq',
#                     height=node.aln_height, width=node.aln_width, scale_factor=node.scale_factor)
#                 ete3.faces.add_face_to_node(SequenceFace, node, 1, aligned=True)
#     else:
#         node.img_style["size"] = 0


# def _get_phyloalignment_colors(root=False):
#         bg = '#000000'
#         fg = '#FFFFFF'
#         bg_colors = {c: bg for c in string.ascii_uppercase}
#         bg_colors['.'] = '#FFFFFF'
#         bg_colors['-'] = '#d3d3d3'
#         fg_colors = {c: fg for c in string.ascii_uppercase}
#         fg_colors['.'] = '#000000'
#         fg_colors['-'] = '#000000'
#         return bg_colors, fg_colors


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
