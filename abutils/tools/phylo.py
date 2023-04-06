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
from typing import Callable, Iterable, Optional, Union
import uuid

import matplotlib.pyplot as plt
from matplotlib import markers
from matplotlib.path import Path

import baltic as bt

from Bio import Phylo

from abstar.core.germline import get_imgt_germlines

from ..core.sequence import Sequence, to_fasta
from ..utils.alignment import mafft
from ..tools.cluster import cluster
from ..utils.pipeline import make_dir


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
        fasttree_bin = os.path.join(
            mod_dir, f"bin/fasttree_{platform.system().lower()}"
        )
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
        ft_cmd = "fasttree {} > {}".format(alignment_file, tree_file)
    else:
        ft_cmd = "fasttree -nt {} > {}".format(alignment_file, tree_file)
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
        self.sequences = sequences
        self.name = name if name is not None else uuid.uuid4()
        self.id_key = (id_key,)
        self.seq_key = (sequence_key,)
        self._root = root
        self._rename = rename
        self.do_clustering = cluster
        self.clustering_threshold = clustering_threshold
        self._fasta_string = None
        self._aln_string = None
        self._tree_string = None
        self._germ_db = None
        self._tree = None
        self._sizes = None
        self._clusters = None

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
        """ """
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
            else:
                sequences = self.sequences
            if self.root is not None:
                sequences = [self.root] + sequences
            self._fasta_string = to_fasta(
                sequences, id_key=self.id_key, sequence_key=self.seq_key, as_string=True
            )
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
        if self._germ_db is None:
            db = Counter(
                [s.get("germline_database", None) for s in self.sequences]
            ).most_common()[0][0]
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
        self._fasta_string = None
        self._aln_string = None
        self._tree_string = None

    def rename(self, old: str) -> str:
        """
        Sequence renamer
        """
        if self._rename is None:
            return old
        if callable(self._rename):
            new = self._rename(old)
        else:
            new = self._rename.get(old, None)
        return new if new is not None else old

    def get_size(self, name, min_size=1, default=1, multiplier=1):
        """ """
        s = self.sizes.get(name, default)
        if s >= min_size:
            return s * multiplier
        return 0

    def cluster(self):
        """
        Cluster sequences prior to alignment
        """
        clusters = cluster(
            self.sequences, threshold=self.clustering_threshold, algo="vsearch", iddef=1
        )
        sizes = {}
        centroids = []
        for c in clusters:
            centroids.append(c.centroid)
            sizes[c.centroid.id] = c.size
        self._sizes = sizes
        self._clusters = clusters

    def plot(
        self,
        figsize=[8, 8],
        size=None,
        color=None,
        min_size=1,
        size_multiplier=10,
        **kwargs,
    ):
        """
        Plot the phylogenetic tree.
        """
        plt.figure(figsize=figsize)
        ax = plt.gca()

        if size is None:
            size = lambda k: self.get_size(
                k.name, default=1, min_size=min_size, multiplier=size_multiplier
            )
        if color is None:
            color = lambda k: "k"
        self.tree.plotTree(ax, x_attr=lambda k: k.height, width=2)
        self.tree.plotPoints(
            ax,
            x_attr=lambda k: k.height,
            size=size,
            colour=color,
            zorder=100,
            outline=0,
            marker=align_marker("o", halign="left"),
            **kwargs,
        )

        # remove spines, ticks and ticklabels
        for s in ax.spines.keys():
            ax.spines[s].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(length=0, labelsize=0)

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
    sequences: Iterable[Sequence],
    name: Optional[str] = None,
    root: Optional[Union[str, Sequence]] = None,
    cluster: bool = True,
    clustering_threshold: float = 1.0,
    rename: Optional[Union[dict, Callable]] = None,
    id_key: Optional[str] = None,
    sequence_key: Optional[str] = None,
) -> Phylogeny:
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
    return Phylogeny(
        sequences=sequences,
        name=name,
        root=root,
        cluster=cluster,
        clustering_threshold=clustering_threshold,
        rename=rename,
        id_key=id_key,
        sequence_key=sequence_key,
    )


def align_marker(
    marker: str = "o", halign: str = "center", valign: str = "middle"
) -> Path:
    """
    create markers with specified alignment.

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
    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = m.get_path().transformed(m.get_transform()).vertices
    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, m.get_path().codes)
