#!/usr/bin/env python
# filename: phylo.py


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


from collections import Counter
from copy import copy, deepcopy
from io import StringIO
import os
import platform
import subprocess as sp
import tempfile
from typing import Callable, Iterable, List, Optional, Union
import uuid

import matplotlib.pyplot as plt
from matplotlib import markers
from matplotlib.path import Path

import baltic as bt

from Bio import Phylo

from abstar.core.germline import get_imgt_germlines

from ..core.sequence import Sequence, to_fasta, read_fasta
from ..utils.alignment import mafft
from ..tools.cluster import cluster
from ..utils.pipeline import make_dir


__all__ = ["fasttree", "Phylogeny", "phylogeny"]


def fasttree(
    aln: str,
    tree_file: Optional[str] = None,
    is_aa: bool = False,
    fasttree_bin: Optional[str] = None,
    debug: bool = False,
    quiet: bool = True,
) -> str:
    """
    Computes a tree file from a multiple seqeunce alignment using `FastTree`_.

    Parameters
    ----------
    aln : str
        Path to a multiple sequence alignment file, in FASTA format, or a
        FASTA-formatted multiple sequence alignment string. Required.

    tree_file : str
        Path to the tree file which will be output by FastTree. If the parent
        directory does not exist, it will be created. If not provided, the output
        (a Newick-formatted tree file) will be returned as a ``str``.

    is_aa : bool, default=False
        Must be set to ``True`` if the input multiple sequence alignment contains
        amino acid sequences. Default is ``False``, meaning FastTree will expect
        nucleotide sequences.

    fasttree_bin : str, optional
        Path to the desired FastTree binary. Default is to use the version of
        FastTree that is bundled with ``abutils``.

    debug : bool, default=False
        If ``True``, verbose output is printed. Default is False.

    quiet : bool, default=True
        Depricated, but retained for backwards compatibility. Use `debug` instead.


    Returns
    -------
    tree_file: str
        Path to the tree file produced by FastTree.


    .. _FastTree:
        http://www.microbesonline.org/fasttree/
    """
    # process input
    if os.path.isfile(aln):
        alignment_file = os.path.abspath(aln)
    else:
        ff = tempfile.NamedTemporaryFile(delete=False)
        ff.close()
        alignment_file = ff.name
        with open(alignment_file, "w") as f:
            f.write(aln)
    if not quiet:
        debug = True
    # set the FastTree binary
    if fasttree_bin is None:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # fasttree_bin = os.path.join(
        #     mod_dir, f"bin/fasttree_{platform.system().lower()}"
        # )
        system = platform.system().lower()
        machine = platform.machine().lower()
        fasttree_bin = os.path.join(mod_dir, f"bin/fasttree_{system}_{machine}")
        fasttree_bin = fasttree_bin.replace("x86_64", "amd64")
    # make output directory if necessary
    if tree_file is None:
        as_file = False
        tree_file = tempfile.NamedTemporaryFile(delete=False).name
    else:
        as_file = True
        tree_file = os.path.abspath(tree_file)
    if not os.path.isdir(os.path.dirname(tree_file)):
        make_dir(os.path.dirname(tree_file))
    # run FastTree
    if is_aa:
        ft_cmd = f"{fasttree_bin} {alignment_file} > {tree_file}"
    else:
        ft_cmd = f"{fasttree_bin} -nt {alignment_file} > {tree_file}"
    ft = sp.Popen(ft_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = ft.communicate()
    if debug:
        print(ft_cmd)
        print(stdout)
        print(stderr)
    # output
    if as_file:
        return tree_file
    else:
        with open(tree_file, "r") as f:
            tree_string = f.read()
        return tree_string


class Phylogeny:
    """
    Base phylogeny class
    """

    def __init__(
        self,
        sequences: Iterable[Sequence],
        name: Optional[str] = None,
        root: Optional[Union[str, Sequence]] = None,
        cluster: bool = True,
        clustering_threshold: float = 1.0,
        clustering_algo: str = "auto",
        rename: Optional[Union[dict, Callable]] = None,
        id_key: Optional[str] = None,
        sequence_key: Optional[str] = None,
    ):
        """
        Phylogenetic representation of an antibody lineage.

        Parameters
        ----------
        sequences : list of Sequence
            A list of ``abutils.Sequence`` objects. Required.

        name : str, default=None
            Name of the lineage. If not provided, a random name will be generated
            using ``uuid.uuid4()``.

        root : str or Sequence, default=None
            Root of the phylogenetic tree. Can be either a sequence ID (the root sequence must
            be in `sequences`) or an ``abutils.Sequence``. If the provided ``Sequence`` is
            already in `sequences`, the duplicate will be ignored when constructing the tree.
            If not provided, the germline V-gene will be used as root.

        cluster : bool, default=True
            Whether or not to cluster seqeunces by identity prior to alignment and tree inference.

        clustering_threshold : float, default=1.0
            Identity threshold for clustering. Must be between 0-1. Default is ``1.0``, which
            collapses identical sequences.

        rename : dict or Callable, default=None
            Used to rename `sequences`. Can be either a ``dict`` of the format:
                ``{old_name: new_name, ...}``
            or a callable function that accepts the old name and returns the new name. Names
            not found in the ``dict`` or for which the function returns ``None`` will not be
            renamed. If not provided, sequences are not renamed.

        id_key : str, default=None
            Key to retrieve the sequence ID. If not provided or missing, ``Sequence.id`` is used.

        sequence_key : str, default=None
            Key to retrieve the sequence. If not provided or missing, ``Sequence.sequence`` is used.

        """
        self.id_key = id_key
        self.seq_key = sequence_key
        self.name = name if name is not None else str(uuid.uuid4())
        self.do_clustering = cluster
        self.clustering_threshold = clustering_threshold
        self.clustering_algo = clustering_algo
        self._root = root
        self._rename = rename
        self._fasta_string = None
        self._aln_string = None
        self._tree_string = None
        self._germ_db = None
        self._tree = None
        self._sizes = None
        self._clusters = None
        self._cluster_dict = None
        # wait until all other attributes are set before processing sequences
        self.sequences = self.process_sequences(sequences)

    def __eq__(self, other):
        """
        ``Phylogeny`` class instances are equal if their names are the same.
        """
        return self.name == other.name

    def __iter__(self):
        """
        Iterates through ``self.tree.Objects``.
        """
        if self._tree is None:
            return None
        for o in self.tree.Objects:
            yield o

    @property
    def tree(self):
        """
        Baltic ``Tree`` object
        """
        if self._tree is None:
            tree_file = StringIO(self.tree_string)
            self._tree = bt.loadNewick(tree_file, tip_regex="_([0-9\-]+)$")
        return self._tree

    @property
    def root(self):
        """
        Root of the tree.
        """
        if self._root is None:
            top_v = Counter(
                [s.get("v_call", None) for s in self.sequences]
            ).most_common()[0][0]
            if top_v is None:
                return
            germ = get_imgt_germlines(self.germ_db, "V", gene=top_v)
            self._root = self._get_top_germline_v()
        elif isinstance(self._root, str):
            root = [s for s in self.sequences if s.id == self._root]
            if not root:
                err = f"WARNING: provided root name ({self._root}) was not found in "
                err += "the list of input sequences. Using the germline V-gene as root."
                print(err)
                self._root = self._get_top_germline_v()
            else:
                self._root = root[0]
        elif self._root is False:
            return None
        return self._root

    @root.setter
    def root(self, root):
        if isinstance(root, str):
            r = [s for s in self.sequences if s.id == root]
            if not r:
                err = f"WARNING: provided root name ({r}) was not found in "
                err += "the list of input sequences. Root was not updated."
                print(err)
            else:
                self._root = r[0]
        else:
            self._root = Sequence(root)

    @property
    def fasta_string(self):
        if self._fasta_string is None:
            if self.do_clustering:
                sequences = self.clusters.centroids
            if self.root is not None:
                sequences = [self.root] + sequences
            self._fasta_string = to_fasta(sequences)
        return self._fasta_string

    @fasta_string.setter
    def fasta_string(self, fasta_string):
        self._fasta_string = fasta_string
        self._clusters = None
        self._aln_string = None
        self._tree_string = None

    @property
    def aln_string(self):
        """
        Multiple sequence alignment, as a FASTA-formatted string.
        """
        if self._aln_string is None:
            aln = mafft(self.fasta_string, as_string=True)
            self._aln_string = aln
        return self._aln_string

    @aln_string.setter
    def aln_string(self, aln_string):
        self._aln_string = aln_string
        self._tree_string = None

    @property
    def tree_string(self):
        """
        Newick tree file, as a string.
        """
        if self._tree_string is None:
            tree_string = fasttree(self.aln_string)
            # reroot only if the root sequence exists
            if self.root is not None:
                tree = Phylo.read(StringIO(tree_string), "newick")
                tree.root_with_outgroup(self.root.id)
                tree_string = tree.format("newick")
            self._tree_string = tree_string
        return self._tree_string

    @tree_string.setter
    def tree_string(self, tree_string):
        self._tree_string = tree_string

    @property
    def germ_db(self):
        """
        Identifies the germline database used to annotate the sequences. If
        multiple databases were used, the most common one is selected.
        """
        if self._germ_db is None:
            db_counts = Counter(
                [s.get("germline_database", None) for s in self.sequences]
            )
            db = db_counts.most_common()[0][0]
            self._germ_db = db
        return self._germ_db

    @germ_db.setter
    def germ_db(self, germ_db):
        self._germ_db = germ_db

    @property
    def sizes(self):
        """
        Returns a ``dict`` of tip sizes using clustering results. Only clusters
        with more than one sequence will be in the ``dict``.
        """
        if self._sizes is None:
            self._sizes = {}
        return self._sizes

    @property
    def clusters(self):
        if self._clusters is None:
            self.cluster()
        return self._clusters

    @clusters.setter
    def clusters(self, clusters):
        self._clusters = clusters
        self._cluster_dict = None
        self._fasta_string = None
        self._aln_string = None
        self._tree_string = None

    @property
    def cluster_dict(self):
        """
        Returns a dict mapping cluster centroid IDs to a list of
        all sequence IDs in the corresponding cluster.
        """
        if self._cluster_dict is None:
            d = {c.centroid.id: c.seq_ids for c in self.clusters}
            self._cluster_dict = d
        return self._cluster_dict

    def sanitize_name(self, name: str) -> str:
        """
        Sanitizes sequence names to be used in the tree, by replacing characters
        that are not allowed in Newick tree files with tokens that can be used to
        un-sanitize the names later.

        Characters that are not allowed in Newick tree files are:
            - colon (:)
            - semicolon (;)
            - comma (,)
            - space ( )
            - open parenthesis (()
            - close parenthesis ())
            - open bracket ([)
            - close bracket (])
            - single quote (')
        """
        return (
            name.replace(":", "<colon>")
            .replace(";", "<semi>")
            .replace(",", "<comma>")
            .replace(" ", "<blank>")
            .replace("(", "<open>")
            .replace(")", "<close>")
            .replace("[", "<openbracket>")
            .replace("]", "<closebracket>")
            .replace("'", "<quote>")
        )

    def unsanitize_name(self, name: str) -> str:
        """
        Un-sanitizes sequence names to be used in the tree, by replacing tokens
        with the original characters.
        """
        return (
            name.replace("<colon>", ":")
            .replace("<semi>", ";")
            .replace("<comma>", ",")
            .replace("<blank>", " ")
            .replace("<open>", "(")
            .replace("<close>", ")")
            .replace("<openbracket>", "[")
            .replace("<closebracket>", "]")
            .replace("<quote>", "'")
        )

    def process_sequences(self, sequences: List[Sequence]) -> List[Sequence]:
        """
        Processes sequences to be used in the tree. Sequence names are sanitized.

        Parameters
        ----------
        sequences : list of Sequence objects
            Sequences to be processed.

        Returns
        -------
        sequences : list of Sequence objects
            Processed sequences.
        """
        sequences = deepcopy(sequences)
        sequences = [
            Sequence(
                s[self.seq_key] if self.seq_key is not None else s.sequence,
                id=self.sanitize_name(
                    self.rename(s[self.id_key] if self.id_key is not None else s.id)
                ),
            )
            for s in sequences
        ]
        return sequences

    def rename(self, old: str) -> str:
        """
        Sequence renamer. The ``old`` name is sanitized before being passed to
        the ``rename`` function. If ``rename`` is a ``dict``, the unsanitized name
        is used as the key to look up the new name. If ``rename`` is a ``callable``,
        the unsanitized name is passed to the function and the return value is used
        as the new name. If ``rename`` is ``None``, the unsanitized name is returned.
        """
        old = self.unsanitize_name(old)
        if self._rename is None:
            return old
        if callable(self._rename):
            new = self._rename(old)
        else:
            new = self._rename.get(old, None)
        return new if new is not None else old

    def get_size(
        self,
        name: str,
        size: Union[Callable, dict, int, float, None] = None,
        min_size: Union[int, float] = 1,
        default: Union[int, float] = 1,
        multiplier: Union[int, float] = 1,
    ) -> Union[int, float]:
        """
        Get the size of a leaf marker, based on the provided ``size`` argument.

        Parameters
        ----------
        name : str
            Name of the leaf node.

        size : callable, dict, int, float, or None
            If ``callable``, it should take a single argument (the name of the
            leaf node) and return a size. If ``dict``, it should map leaf node
            names to sizes. If ``int`` or ``float``, it should be the size of
            all leaf nodes. If ``None``, the size of all leaf nodes will be 1.

        min_size : int or float
            Minimum size of a leaf node. If the size of a leaf node is less than
            ``min_size``, it will be set to 0.

        default : int or float
            Default size of a leaf node, if it is not found in the ``size``
            argument.

        multiplier : int or float
            Multiplier for the size of a leaf node.

        Returns
        -------
        size : int or float
            Size of the leaf node.
        """
        unsanitized_name = self.unsanitize_name(name)
        if callable(size):
            s = size(unsanitized_name)
        elif isinstance(size, dict):
            s = size.get(unsanitized_name, default)
        elif isinstance(size, (int, float)):
            s = size
        else:
            s = self.sizes.get(name, default)
        if s >= min_size:
            return s * multiplier
        return 0

    def get_color(
        self,
        name: str,
        color: Union[Callable, dict, str, Iterable, None] = None,
        default: Union[str, Iterable] = "black",
    ) -> Union[str, Iterable]:
        """
        Get the color of a leaf marker, based on the provided ``color`` argument.
        If the provided name is a cluster centroid, all sequences in the corresponding
        cluster will be evaluated and the most common color will be used.

        Parameters
        ----------
        name : str
            Name of the leaf node.

        color : callable, dict, str, iterable, or None
            If ``callable``, it should take a single argument (the name of the
            leaf node) and return a color. If ``dict``, it should map leaf node
            names to colors. If ``str``, it should be the color of all leaf nodes.
            If ``iterable``, it should be the color of all leaf nodes, as an
            iterable of RGB(A) values. If ``None``, the color of all leaf nodes
            will be `default`.

        default : str or iterable
            Default color of a leaf node, if it is not found in the ``color``.

        Returns
        -------
        color : str or iterable
            The color of the marker, as a string or an iterable of RGB(A) values.
        """
        names = self.cluster_dict.get(name, [name])
        colors = []
        for n in names:
            n = self.unsanitize_name(n)
            # color is a function
            if callable(color):
                c = color(n)
            # color is a dict
            elif isinstance(color, dict):
                c = color.get(n, default)
            # a single color, as a string
            elif isinstance(color, str):
                c = color
            # a single color, as an iterable of RGB(A) values
            elif (
                isinstance(color, (list, tuple))
                and 3 <= len(color) <= 4
                and all([isinstance(c, (int, float)) for c in color])
            ):
                c = color
            else:
                c = default
            colors.append(c)
        color = Counter(colors).most_common()[0][0]
        return color

    # TODO: add support for branch color orders
    # useful for things like longitudinal data, where the timepoint of the
    # branch can be used to color by the earliest timepoint in the branch

    def get_branch_color(
        self,
        k: Union[bt.node, bt.leaf],
        color: Union[Callable, dict, str, Iterable, None] = None,
        default: Union[str, Iterable] = "black",
    ) -> Union[str, Iterable]:
        """
        Get the color of a branch, based on the provided ``color`` argument.

        Parameters
        ----------
        k : bt.node or bt.leaf
            The branch to get the color of. If a node, the most common color of
            all child leaves will be returned.

        color : callable, dict, str, iterable, or None
            If ``callable``, it should take a single argument (the name of the
            leaf node) and return a color. If ``dict``, it should map leaf node
            names to colors. If ``str``, it should be the color of all leaf nodes.
            If ``iterable``, it should be the color of all leaf nodes, as an
            iterable of RGB(A) values. If ``None``, the color of all leaf nodes
            will be `default`.

        default : str or iterable
            Default color of a leaf node, if it is not found in the ``color``.

        Returns
        -------
        color : str or iterable
            The color of the branch, as a string or an iterable of RGB(A) values.
        """
        if k.is_node():
            leaves = list(k.leaves)
        else:
            leaves = [k]
        colors = []
        for leaf in leaves:
            if not isinstance(leaf, str):
                leaf = leaf.name
            names = self.cluster_dict.get(leaf, [leaf])
            for name in names:
                c = self.get_color(name, color=color, default=default)
                colors.append(c)
        color = Counter(colors).most_common()[0][0]
        return color

    def get_marker_edgewidth(
        self,
        name: str,
        edgewidth: Union[Callable, dict, int, float, None] = None,
        default: Union[int, float] = 1,
    ) -> Union[int, float]:
        """
        Get the edge width of a leaf marker, based on the provided ``edgewidth`` argument.

        Parameters
        ----------
        name : str
            Name of the leaf node.

        edgewidth : callable, dict, int, float, or None
            If ``callable``, it should take a single argument (the name of the
            leaf node) and return an edge width. If ``dict``, it should map leaf node
            names to edge widths. If ``int`` or ``float``, it should be the edge width of
            all leaf nodes. If ``None``, the edge width of all leaf nodes will be `default`.

        default : int or float
            Default edge width of a leaf node, if it is not found in the ``edgewidth``.

        Returns
        -------
        edgewidth : int or float
            The edge width of the marker.
        """
        name = self.unsanitize_name(name)
        if callable(edgewidth):
            w = edgewidth(name)
        elif isinstance(edgewidth, dict):
            w = edgewidth.get(name, default)
        elif isinstance(edgewidth, (int, float)):
            w = edgewidth
        else:
            w = default
        return w

    def cluster(self):
        """
        Cluster sequences prior to alignment.
        """
        clusters = cluster(
            self.sequences,
            threshold=self.clustering_threshold,
            algo=self.clustering_algo,
            iddef=1,
            id_key=self.id_key,
            seq_key=self.seq_key,
        )
        sizes = {}
        for c in clusters:
            sizes[c.centroid.id] = c.size
        self._clusters = clusters
        self._sizes = sizes

    def plot(
        self,
        size: Union[Callable, dict, int, float, None] = None,
        color: Union[Callable, dict, Iterable, str, None] = None,
        alpha: float = 0.75,
        min_size: int = 1,
        size_multiplier: Union[int, float] = 10,
        linewidth: Union[int, float] = 2,
        color_branches: bool = False,
        x_attr: Callable = lambda k: k.height,
        y_attr: Callable = lambda k: k.y,
        connection_type: str = "baltic",
        radial: bool = False,
        radial_start: float = 0,
        radial_fraction: float = 1.0,
        inward_space: float = 0.1,
        marker: str = "o",
        marker_edgewidth: Union[Callable, dict, int, float] = 0,
        marker_edgecolor: Union[Callable, dict, str, Iterable, None] = None,
        marker_halign: str = "left",
        marker_valign: str = "center",
        figsize: Iterable[Union[int, float]] = [8, 8],
        show: bool = False,
        figfile: Union[str, Path] = None,
        **kwargs,
    ):
        """
        Plot the phylogenetic tree.

        Parameters
        ----------
        size : callable or dict or int or float or None, optional
            Size of the markers at leaf edges. If a callable, it should take
            a ``baltic`` leaf as input and return a size. If a ``dict``, it
            should have leaf names as keys and sizes as values. If an ``int``,
            or ``float``, all markers will have the same size. If ``None``, the
            cluster sizes will be used.

        color : callable or dict or iterable or str or None, optional
            Color of the markers at leaf edges. If a callable, it should take
            a ``baltic`` leaf object as input and return a color. If a ``dict``,
            it should have leaf names as keys and colors as values. If a ``str``,
            or an iterable of RGB(A) values, all markers will have the same color.
            If ``None``, all markers will be black.

        alpha : float, optional
            Alpha value for the markers.

        min_size : int, optional
            Minimum size of the markers. Any leaf edges with a size smaller than
            `min_size` will not have a marker.

        size_multiplier : int or float, optional
            Multiplier for the marker sizes.

        linewidth : int or float, optional
            Width of the tree lines.

        color_branches : bool, optional
            Whether to color the branches of the tree. If ``True``, the branches
            will be colored according to the color of the leaf edges.

        x_attr : callable, optional
            Attribute of the ``baltic`` tree object to use for the x-axis.

        y_attr : callable, optional
            Attribute of the ``baltic`` tree object to use for the y-axis.

        connection_type : str, optional
            Type of connection to use. One of ``'baltic'``, ``'direct'``, or ``'elbow'``.

        radial : bool, optional
            Whether to plot the tree in a radial fashion.

        radial_start : float, optional
            Starting point for the radial plot, as a fraction of the circle.

        radial_fraction : float, optional
            Fraction of the circle to use for the radial plot.

        inward_space : float, optional
            Fraction of the circle to leave empty at the center of the radial plot.

        marker : str, optional
            Marker style. See ``matplotlib.pyplot.scatter`` for more details.

         marker_edgewidth : callable or dict or int or float, optional
            Width of the marker edges. If a callable, it should take a ``baltic`` leaf
            as input and return a width. If a ``dict``, it should have leaf names as
            keys and widths as values. If an ``int``, or ``float``, all markers will
            have the same width.

        marker_edgecolor : callable or dict or str or iterable or None, optional
            Color of the marker edges. If a callable, it should take a ``baltic`` leaf
            object as input and return a color. If a ``dict``, it should have leaf names
            as keys and colors as values. If a ``str``, or an iterable of RGB(A) values,
            all markers will have the same color. If ``None``, all marker edges will
            match the marker color.

        marker_halign : str, optional
            Horizontal alignment of the markers. One of ``'left'``, ``'center'``, or
            ``'right'``.

        marker_valign : str, optional
            Vertical alignment of the markers. One of ``'top'``, ``'center'``, or
            ``'bottom'``.

        figsize : iterable of int or float, optional
            Figure size, by default [8, 8]

        show : bool, optional
            Whether to show the figure.

        figfile : str or Path, optional
            Path to save the figure to.

        """
        plt.figure(figsize=figsize)
        ax = plt.gca()
        # plot parameters
        size_func = lambda k: self.get_size(
            k.name,
            size=size,
            default=1,
            min_size=min_size,
            multiplier=size_multiplier,
        )
        color_func = lambda k: self.get_color(k.name, color=color, default="black")
        branch_color_func = lambda k: self.get_branch_color(
            k, color=color, default="black"
        )
        marker = align_marker(marker, halign=marker_halign, valign=marker_valign)
        if marker_edgecolor is None:
            marker_edgecolor_func = color_func
        else:
            marker_edgecolor_func = lambda k: self.get_color(
                k.name, color=marker_edgecolor, default="black"
            )
        marker_edgewidth_func = lambda k: self.get_marker_edgewidth(
            k.name, edgewidth=marker_edgewidth, default=0
        )
        # make the plot
        if radial:
            # circular phlogenetic tree
            inward_space = inward_space * self.tree.treeHeight
            self.tree.plotCircularTree(
                ax,
                x_attr=x_attr,
                y_attr=y_attr,
                width=linewidth,
                colour=branch_color_func if color_branches else None,
                connection_type=connection_type,
                circStart=radial_start,
                circFracn=radial_fraction,
                inwardSpace=inward_space,
            )
            self.tree.plotCircularPoints(
                ax,
                x_attr=x_attr,
                y_attr=y_attr,
                size=size_func,
                colour=color_func,
                zorder=100,
                marker=marker,
                outline_size=marker_edgewidth_func,
                outline_colour=marker_edgecolor_func,
                circStart=radial_start,
                circFrac=radial_fraction,
                inwardSpace=inward_space,
                alpha=alpha,
                **kwargs,
            )
        else:
            # standard phlogenetic tree
            self.tree.plotTree(
                ax,
                x_attr=x_attr,
                y_attr=y_attr,
                width=linewidth,
                colour=branch_color_func if color_branches else None,
                connection_type=connection_type,
            )
            self.tree.plotPoints(
                ax,
                x_attr=x_attr,
                y_attr=y_attr,
                size=size_func,
                colour=color_func,
                zorder=100,
                marker=marker,
                outline_size=marker_edgewidth_func,
                outline_colour=marker_edgecolor_func,
                alpha=alpha,
                **kwargs,
            )
        # style the plot
        for s in ax.spines.keys():
            ax.spines[s].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(length=0, labelsize=0)
        # save or show
        if show:
            plt.tight_layout()
            plt.show()
        elif figfile is not None:
            plt.tight_layout()
            plt.savefig(figfile)
        else:
            return ax

    def _get_top_germline_v(self) -> Union[Sequence, None]:
        """
        Retrieves the germline V-gene sequence. If sequences were assigned different V-gene
        alleles, the most common V-gene allele is used.

        The appropriate germline database is inferred from sequence annotations.

        If the germline V-gene can't be found, the germline database can't be found/inferred,
        or the sequences don't have AIRR-compatible annotations, returns ``None``.
        """
        v_counts = Counter([s.get("v_call", None) for s in self.sequences])
        top_v = v_counts.most_common()[0][0]
        if top_v is None:
            return None
        try:
            germ = get_imgt_germlines(self.germ_db, "V", gene=top_v)
            return Sequence(germ.ungapped_nt_sequence, id=top_v)
        except:
            return


def phylogeny(
    sequences: Union[str, Iterable[Sequence]],
    name: Optional[str] = None,
    root: Optional[Union[str, Sequence]] = None,
    cluster: bool = True,
    clustering_threshold: float = 1.0,
    clustering_algo: str = "auto",
    rename: Optional[Union[dict, Callable]] = None,
    id_key: Optional[str] = None,
    sequence_key: Optional[str] = None,
) -> Phylogeny:
    """
    Phylogenetic representation of an antibody lineage.

    Parameters
    ----------
    sequences : str or list of Sequence
        A list of ``abutils.Sequence`` objects or the path to a FASTA-formatted file.
        Required.

    name : str, default=None
        Name of the lineage. If not provided, a random name will be generated
        using ``uuid.uuid4()``.

    root : str or Sequence, default=None
        Root of the phylogenetic tree. Can be either a sequence ID (the root sequence must
        be in `sequences`) or an ``abutils.Sequence``. If the provided ``Sequence`` is
        already in `sequences`, the duplicate will be ignored when constructing the tree.
        If not provided, the germline V-gene will be used as root.

    cluster : bool, default=True
        Whether or not to cluster seqeunces by identity prior to alignment and tree inference.

    clustering_threshold : float, default=1.0
        Identity threshold for clustering. Must be between 0-1. Default is ``1.0``, which
        collapses identical sequences.

    rename : dict or Callable, default=None
        Used to rename `sequences`. Can be either a ``dict`` of the format:
            ``{old_name: new_name, ...}``
        or a callable function that accepts the old name and returns the new name. Names
        not found in the ``dict`` or for which the function returns ``None`` will not be
        renamed. If not provided, sequences are not renamed.

    id_key : str, default=None
        Key to retrieve the sequence ID. If not provided or missing, ``Sequence.id`` is used.

    sequence_key : str, default=None
        Key to retrieve the sequence. If not provided or missing, ``Sequence.sequence`` is used.

    Returns
    -------
    phylogeny : Phylogeny
        An ``abutils.Phylogeny`` object.
    """
    if isinstance(sequences, str):
        sequences = read_fasta(sequences)
    return Phylogeny(
        sequences=sequences,
        name=name,
        root=root,
        cluster=cluster,
        clustering_threshold=clustering_threshold,
        clustering_algo=clustering_algo,
        rename=rename,
        id_key=id_key,
        sequence_key=sequence_key,
    )


def align_marker(
    marker: str = "o", halign: str = "center", valign: str = "middle"
) -> Path:
    """
    Create markers with specified alignment.
    https://stackoverflow.com/questions/26686722/align-matplotlib-scatter-marker-left-and-or-right


    Parameters
    ----------
    marker : a valid marker specification.
      See mpl.markers

    halign : string, float {'left', 'center', 'right'}
      Specifies the horizontal alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'center',
      -1 is 'right', 1 is 'left').

    valign : string, float {'top', 'middle', 'bottom'}
      Specifies the vertical alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'middle',
      -1 is 'top', 1 is 'bottom').

    Returns
    -------
    marker_array : numpy.ndarray
      A Nx2 array that specifies the marker path relative to the
      plot target point at (0, 0).

    Notes
    -----
    The mark_array can be passed directly to ax.plot and ax.scatter, e.g.::

        ax.plot(1, 1, marker=align_marker('>', 'left'))

    """
    if isinstance(halign, str):
        halign = {
            "right": -1.0,
            "middle": 0.0,
            "center": 0.0,
            "left": 1.0,
        }[halign]
    if isinstance(valign, str):
        valign = {
            "top": -1.0,
            "middle": 0.0,
            "center": 0.0,
            "bottom": 1.0,
        }[valign]
    m = markers.MarkerStyle(marker)
    m_arr = m.get_path().transformed(m.get_transform()).vertices
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2
    return Path(m_arr, m.get_path().codes)
