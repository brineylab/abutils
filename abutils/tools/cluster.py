#!/usr/bin/env python
# filename: cluster.py


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

from __future__ import annotations

import os

# import platform
import random
import string
import subprocess as sp

# import sys
import tempfile
from typing import Iterable

# from Bio.Align import AlignInfo
from ..bin import get_path as get_binary_path
from ..core.pair import Pair
from ..core.sequence import Sequence
from ..io import read_fasta, to_fasta

# from ..utils.decorators import lazy_property
# from ..utils.pipeline import make_dir
from .alignment import make_consensus

__all__ = [
    "cluster",
    "cluster_vsearch",
    "cluster_mmseqs",
    "cluster_cdhit",
    "linclust",
    "create_mmseqs_db",
    "Clusters",
    "Cluster",
]


# -----------------
# Helper functions
# -----------------


def _is_mmseqs_db(path: str) -> bool:
    """
    Check if path points to an existing MMseqs2 database.

    MMseqs2 databases consist of multiple files with a common prefix.
    The main database file has no extension, and auxiliary files
    include .index, .dbtype, etc.

    Parameters
    ----------
    path : str
        Path to check.

    Returns
    -------
    bool
        True if the path points to a valid MMseqs2 database.
    """
    if not isinstance(path, str):
        return False
    # MMseqs2 databases have .index and .dbtype files
    required_files = [path, f"{path}.index", f"{path}.dbtype"]
    return all(os.path.exists(f) for f in required_files)


def _validate_output_path(path: str) -> None:
    """
    Validate that output path's parent directory exists.

    Parameters
    ----------
    path : str
        Output file path to validate.

    Raises
    ------
    FileNotFoundError
        If the parent directory does not exist.
    """
    parent_dir = os.path.dirname(os.path.abspath(path))
    if parent_dir and not os.path.exists(parent_dir):
        raise FileNotFoundError(f"Parent directory does not exist: {parent_dir}")


def _cleanup_mmseqs_db(db_path: str) -> None:
    """
    Remove all files associated with an MMseqs2 database.

    Parameters
    ----------
    db_path : str
        Path prefix of the MMseqs2 database to clean up.
    """
    extensions = [
        "",
        ".index",
        ".dbtype",
        "_h",
        "_h.index",
        "_h.dbtype",
        ".lookup",
        ".source",
    ]
    for ext in extensions:
        f = f"{db_path}{ext}"
        if os.path.exists(f):
            os.remove(f)


class Cluster:
    """A sequence cluster from clustering operations.

    Represents a group of similar sequences identified through clustering.
    Provides access to cluster members, centroid, and consensus sequences.

    Args:
        name: Cluster identifier.
        sequences: List of Sequence objects in the cluster.
        centroid: Centroid Sequence of the cluster. Defaults to ``None``.

    Attributes:
        name: Cluster identifier.
        sequences: List of Sequence objects in the cluster.
        size: Number of sequences in the cluster.
        sequence_ids: List of sequence IDs in the cluster.
        centroid: Centroid Sequence (most representative member).
        consensus: Consensus Sequence computed from cluster alignment.
    """

    def __init__(self, name, sequences, centroid=None):
        self.name = name
        self.sequences = sequences
        self.centroid = centroid
        self._consensus = None

    def __iter__(self):
        for s in self.sequences:
            yield s

    @property
    def size(self) -> int:
        return len(self.sequences)

    @property
    def sequence_ids(self) -> Iterable[str]:
        return [s.id for s in self.sequences]

    @property
    def seq_ids(self) -> Iterable[str]:
        # retained for backwards compatibility, use sequence_ids instead
        return [s.id for s in self.sequences]

    @property
    def consensus(self) -> Sequence:
        if self._consensus is None:
            self.make_consensus()
        return self._consensus()

    def make_consensus(
        self,
        threshold: float = 0.51,
        algo: str = "mafft",
        ambiguous: str | None = None,
    ) -> Sequence:
        """Compute a consensus sequence from cluster members.

        Aligns cluster sequences and generates a consensus based on the most
        common base/residue at each position.

        Args:
            threshold: Minimum frequency (0-1) for calling a consensus position.
                Defaults to ``0.51``.
            algo: Multiple sequence alignment algorithm. Options are ``"mafft"``,
                ``"famsa"``, or ``"muscle"``. Defaults to ``"mafft"``.
            ambiguous: Character for ambiguous positions. If ``None``, uses the
                default for the sequence type. Defaults to ``None``.

        Returns:
            Consensus Sequence object.
        """
        if len(self.sequences) == 1:
            self._consensus = self.sequences[0]
        else:
            consensus = make_consensus(
                self.sequences,
                threshold=threshold,
                algo=algo,
                ambiguous=ambiguous,
            )
            self._consensus = consensus
        return self._consensus
        # return consensus
        # aln = mafft(self.sequences)
        # if aln is None:
        #     print("ERROR: Failed to generate an alignmnet for a consensus sequence.")
        #     return None
        # summary_align = AlignInfo.SummaryInfo(aln)
        # consensus = summary_align.gap_consensus(threshold=0.51, ambiguous="n")
        # consensus_string = str(consensus).replace("-", "")
        # consensus_seq = Sequence(consensus_string.upper())
        # return consensus_seq


class Clusters:
    """Collection of sequence clusters.

    Container for multiple Cluster objects, providing iteration, indexing,
    and access to cluster statistics.

    Args:
        clusters: Initial clusters as a list of Cluster objects, a dictionary
            mapping cluster names to cluster data, or ``None`` for an empty
            collection. Defaults to ``None``.

    Attributes:
        clusters: List of Cluster objects, sorted by size (largest first).
        centroids: List of centroid Sequences from all clusters.
        largest_cluster: The Cluster with the most sequences.
        count: Number of clusters in the collection.
    """

    def __init__(self, clusters=None):
        self._clusters = self._parse_clusters(clusters)

    def __iter__(self):
        for cluster in self.clusters:
            yield cluster

    def __getitem__(self, key):
        return self.clusters[key]

    def __len__(self):
        return len(self._clusters)

    @property
    def clusters(self):
        return sorted(self._clusters, key=lambda x: x.size, reverse=True)

    @property
    def centroids(self):
        return [c.centroid for c in self.clusters]

    @property
    def largest_cluster(self):
        return self.clusters[0]

    @property
    def count(self):
        return len(self.clusters)

    def _parse_clusters(self, clusts):
        if isinstance(clusts, dict):
            clusters = []
            for name, cdata in clusts.items():
                seqs = cdata.get("seqs", [])
                centroid = cdata.get("centroid", None)
                clusters.append(Cluster(name, seqs, centroid=centroid))
            return clusters
        elif clusts is None:
            return []
        else:
            return clusts

    def add(self, cluster):
        self._clusters.append(cluster)


def cluster(
    sequences: Iterable | str,
    threshold: float = 0.975,
    algo: str = "auto",
    temp_dir: str = "/tmp",
    iddef: int = 0,
    cluster_mode: str = "2",
    cov_mode: str = "0",
    coverage: float = 0.8,
    alignment_mode: str = "3",
    seq_id_mode: str = "1",
    vsearch_bin: str = None,
    mmseqs_bin: str = None,
    cdhit_bin: str = None,
    id_key: str | None = None,
    seq_key: str | None = None,
    threads: int | None = None,
    strand: str = "plus",
    as_dict: bool = False,
    quiet: bool = False,
    debug: bool = False,
) -> dict | Clusters:
    """Cluster sequences by identity using CD-HIT, VSEARCH, or MMseqs2.

    Groups sequences by similarity using one of several clustering algorithms.
    By default, automatically selects the algorithm based on sequence count
    and type.

    Args:
        sequences: Input sequences as Sequence objects, FASTA string/file path,
            BioPython SeqRecords, or list of [id, sequence] pairs.
        threshold: Identity threshold for clustering (0-1). Defaults to ``0.975``.
        algo: Clustering algorithm. Options: ``"vsearch"``, ``"mmseqs"``,
            ``"cdhit"``, or ``"auto"``. Auto selects VSEARCH for <10k nucleotide
            sequences, CD-HIT for <10k amino acid sequences, and MMseqs2
            otherwise. Defaults to ``"auto"``.
        temp_dir: Directory for temporary clustering files. Defaults to ``"/tmp"``.
        iddef: VSEARCH identity definition (0-4). Defaults to ``0``.
        cluster_mode: MMseqs2 clustering mode. ``"1"`` greedy set cover,
            ``"2"`` connected component, ``"3"`` greedy incremental.
            Defaults to ``"2"``.
        cov_mode: MMseqs2 coverage mode (0-3). Defaults to ``"0"``.
        coverage: MMseqs2 coverage threshold (0-1). Defaults to ``0.8``.
        alignment_mode: MMseqs2 alignment mode (0-4). Defaults to ``"3"``.
        seq_id_mode: MMseqs2 sequence ID mode (0-2). Defaults to ``"1"``.
        vsearch_bin: Path to VSEARCH binary. Uses bundled binary if ``None``.
        mmseqs_bin: Path to MMseqs2 binary. Uses bundled binary if ``None``.
        cdhit_bin: Path to CD-HIT binary. Uses bundled binary if ``None``.
        id_key: Key for sequence ID in annotations. Defaults to ``None``.
        seq_key: Key for sequence in annotations. Defaults to ``None``.
        threads: Number of threads for clustering. Defaults to ``None``.
        strand: Sequence strand, ``"plus"`` or ``"both"``. Defaults to ``"plus"``.
        as_dict: If ``True``, return dict instead of Clusters object.
            Defaults to ``False``.
        quiet: If ``True``, suppress clustering output. Defaults to ``False``.
        debug: If ``True``, print debug information. Defaults to ``False``.

    Returns:
        Clusters object or dict mapping cluster names to centroid and sequences.

    Example:
        >>> seqs = abutils.io.read_fasta("sequences.fasta")
        >>> clusters = abutils.tl.cluster(seqs, threshold=0.97)
        >>> print(f"Found {clusters.count} clusters")
    """
    # check input data to get the number of sequences
    fasta_file = to_fasta(
        sequences, tempfile_dir=temp_dir, id_key=id_key, sequence_key=seq_key
    )
    seqs = read_fasta(fasta_file)
    # check for AA sequences (which can't be handled by USEARCH)
    nt_residues = set(["A", "C", "G", "T", "N", "U", "-"])
    is_aa = False
    for s in seqs[:25]:
        if any([r not in nt_residues for r in s.sequence]):
            is_aa = True
    seq_dict = {s.id: s for s in seqs}
    # select the clustering algo
    algo = algo.lower()
    if algo == "auto" and len(seqs) < 10000:
        if is_aa:
            algo = "cdhit"
        else:
            algo = "vsearch"
    elif algo == "auto" and len(seqs) >= 10000:
        algo = "mmseqs"
    if algo in ["mmseqs", "mmseqs2"]:
        cluster_dict = cluster_mmseqs(
            fasta_file=fasta_file,
            threshold=threshold,
            cluster_mode=cluster_mode,
            cov_mode=cov_mode,
            coverage=coverage,
            alignment_mode=alignment_mode,
            seq_id_mode=seq_id_mode,
            threads=threads,
            temp_dir=temp_dir,
            mmseqs_bin=mmseqs_bin,
            as_dict=True,
            quiet=quiet,
            debug=debug,
        )
    elif algo in ["cdhit", "cd_hit", "cd-hit"]:
        cluster_dict = cluster_cdhit(
            fasta_file=fasta_file,
            threshold=threshold,
            temp_dir=temp_dir,
            cdhit_bin=cdhit_bin,
            threads=0 if threads is None else threads,
            as_dict=True,
            quiet=quiet,
            debug=debug,
        )
    elif algo == "vsearch":
        cluster_dict = cluster_vsearch(
            fasta_file=fasta_file,
            threshold=threshold,
            temp_dir=temp_dir,
            iddef=iddef,
            vsearch_bin=vsearch_bin,
            strand=strand,
            threads=threads,
            as_dict=True,
            quiet=quiet,
            debug=debug,
        )
    else:
        err = f"Invalid algo option: {algo}."
        err += " Valid choices are: 'vsearch', 'mmseqs', 'cdhit', or 'auto'."
        raise (ValueError, err)
    cluster_info = {}
    for cname, cdata in cluster_dict.items():
        cluster_info[cname] = {
            "centroid": seq_dict[cdata["centroid_id"]],
            "seqs": [seq_dict[seq_id] for seq_id in cdata["seq_ids"]],
        }
    if as_dict:
        return cluster_info
    return Clusters(cluster_info)


def cluster_vsearch(
    fasta_file: str,
    threshold: float = 0.975,
    temp_dir: str = "/tmp",
    iddef: int = 0,
    vsearch_bin: str = None,
    strand: str = "plus",
    threads: int | None = None,
    as_dict: bool = False,
    quiet: bool = False,
    debug: bool = False,
) -> dict | Clusters:
    """
    Clusters sequences using `VSEARCH <https://github.com/torognes/vsearch>`_.

    Parameters
    ----------
    fasta_file : string
        Path to a FASTA-formatted file. Required. If you'd like to run ``vsearch``
        using ``Sequence`` objects as input, use ``cluster(algo="vsearch")``.

    threshold : float, default=0.975
        Identity threshold for clustering. Must be between 0 and 1.

    temp_dir : str, default="/tmp"
        Path to a directory for temporary storage of clustering files.

    iddef : int, default=1
        Identity definition, as implemented by VSEARCH. Options are:
            0. CD-HIT definition: (matching columns) / (shortest sequence length).
            1. edit distance: (matching columns) / (alignment length).
            2. edit distance excluding terminal gaps (same as --id).
            3. Marine Biological Lab definition counting each gap opening (internal or
               terminal) as a single mismatch, whether or not the gap was extended: 1.0
               - [(mismatches + gap openings)/(longest sequence length)]
            4. BLAST definition, equivalent to --iddef 1 in a context of global pairwise
               alignment.

    vsearch_bin : str, optional
        Path to a VSEARCH executable. If not provided, the VSEARCH binary bundled
        with ``abutils`` will be used.

    id_key : str, default=None
        Key to retrieve the sequence ID. If not provided or missing, ``Sequence.id`` is used.

    sequence_key : str, default=None
        Key to retrieve the sequence. If not provided or missing, ``Sequence.sequence`` is used.

    strand : str, default="plus"
        Strand of the sequences to align. Options are ``"plus"`` and ``"both"``.

    as_dict : bool, default=False
        If ``True``, return clustering results as a ``dict`` rather than a ``Clusters``
        object. the ``dict`` is of the format:
            {"cluster1_name": {"centroid": cluster1_centroid,
                               "seqs": [seq1, seq2, seq3, ...]},
             "cluster2_name": {"centroid": cluster2_centroid,
                               "seqs": [seq4, seq5, seq6, ...]},
            }

    quiet : bool, default=False
        If ``True``, suppresses all output from the clustering algorithm.

    debug : bool, default=False
        If ``True``, prints standard output and standard error from ``vsearch``.
        Default is ``False``.


    Returns
    -------
    clusters : Path to the UC output file from ``vsearch`` or a ``dict`` of cluster info.

    """
    # output files
    centroid_file = tempfile.NamedTemporaryFile(
        dir=temp_dir, delete=False, prefix="centroids_"
    ).name
    uc_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, prefix="uc_").name
    # get the vsearch binary
    if vsearch_bin is None:
        vsearch_bin = get_binary_path("vsearch")
        # mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # system = platform.system().lower()
        # machine = platform.machine().lower().replace("x86_64", "amd64")
        # vsearch_bin = os.path.join(mod_dir, f"bin/vsearch_{system}_{machine}")
    # do clustering
    vsearch_cmd = f"{vsearch_bin} --cluster_fast {fasta_file}"
    vsearch_cmd += f" --centroids {centroid_file}"
    vsearch_cmd += " --clusterout_id"
    vsearch_cmd += f" --uc {uc_file}"
    vsearch_cmd += f" --id {threshold}"
    vsearch_cmd += f" --iddef {iddef}"
    vsearch_cmd += " --sizeout"
    vsearch_cmd += f" --strand {strand}"
    if threads is not None:
        vsearch_cmd += f" --threads {threads}"
    p = sp.Popen(vsearch_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        print("STDOUT:", stdout.decode("utf-8"))
        print("")
        print("STDERR:", stderr.decode("utf-8"))
    # process output
    if os.stat(uc_file).st_size == 0:
        if not quiet:
            err = f"WARNING: the VSEARCH output file ({uc_file}) is empty. "
            err += "Please verify that the input data is valid."
            print(err)
        return {}
    if as_dict:
        cluster_info = {}
        with open(uc_file, "r") as f:
            for line in f:
                if not (ldata := line.strip().split()):
                    continue
                cluster_num = ldata[1]
                if cluster_num not in cluster_info:
                    cluster_info[cluster_num] = {"seq_ids": []}
                if ldata[0] == "S":
                    centroid = ldata[-2]
                    cluster_info[cluster_num]["centroid_id"] = centroid
                    cluster_info[cluster_num]["seq_ids"].append(centroid)
                elif ldata[0] == "H":
                    hit = ldata[-2]
                    cluster_info[cluster_num]["seq_ids"].append(hit)
        if not debug:
            os.remove(uc_file)
            os.remove(centroid_file)
        return cluster_info
    else:
        return uc_file


def cluster_mmseqs(
    fasta_file: str | None = None,
    db_path: str | None = None,
    threshold: float = 0.975,
    cluster_mode: str = "2",
    cov_mode: str = "0",
    coverage: float = 0.8,
    alignment_mode: str = "3",
    seq_id_mode: str = "1",
    threads: int | None = None,
    temp_dir: str = "/tmp",
    mmseqs_bin: str = None,
    as_dict: bool = False,
    quiet: bool = False,
    debug: bool = False,
):
    """
    Clusters sequences using `MMseqs2 <https://github.com/soedinglab/MMseqs2>`_.

    Parameters
    ----------
    fasta_file : string, optional
        Path to a FASTA-formatted file. Either ``fasta_file`` or ``db_path`` must be
        provided. If you'd like to run ``mmseqs`` using ``Sequence`` objects as input,
        use ``cluster(algo="mmseqs")``.

    db_path : string, optional
        Path to a pre-generated MMseqs2 database (created with ``create_mmseqs_db()``).
        Either ``fasta_file`` or ``db_path`` must be provided. Using a pre-generated
        database avoids redundant database creation when clustering the same sequences
        multiple times with different parameters.

    threshold : float, default=0.975
        Identity threshold for clustering. Must be between 0 and 1.

    cluster_mode : str, default="2"
        Clustering mode. Options are ``"1"``, ``"2"``, or ``"3"``.
          - ``"1"``: greedy set cover
          - ``"2"``: connected component
          - ``"3"``: greedy incremental (CD-HIT-like)
        See the `MMseqs2 documentation <https://mmseqs.com/latest/userguide.html#clustering>`_
        for more information.

    cov_mode : str, default="0"
        Coverage mode. Options are ``"0"``, ``"1"``, ``"2"``, or ``"3"``.
            - ``"0"``: bidirectional
            - ``"1"``: target coverage
            - ``"2"``: query coverage
            - ``"3"``: target-in-query length coverage
        See the `MMseqs2 documentation <https://mmseqs.com/latest/userguide.html#clustering>`_
        for more information.

    coverage : float, default=0.8
        Coverage threshold for clustering. Must be between 0 and 1.

    alignment_mode : str, default="3"
        Alignment mode. Options are ``"0"``, ``"1"``, ``"2"``, ``"3"``, or ``"4"``.
            - ``"0"``: automatic
            - ``"1"``: only score and end_pos
            - ``"2"``: also start_pos and cov
            - ``"3"``: also seq.id
            - ``"4"``: only ungapped alignment
        See the `MMseqs2 documentation <https://mmseqs.com/latest/userguide.html#clustering>`_
        for more information.

    seq_id_mode : str, default="1"
        Sequence ID mode. Options are ``"0"``, ``"1"``, or ``"2"``.
            - ``"0"``: alignment length
            - ``"1"``: shorter sequence length
            - ``"2"``: longer sequence length

    threads : int, default=None
        Number of threads to use for clustering. If not provided, the number of threads
        will be determined by MMseqs2.

    temp_dir : str, default="/tmp"
        Path to a directory for temporary storage of clustering files.

    mmseqs_bin : str, optional
        Path to a MMseqs2 executable. If not provided, the MMseqs2 binary bundled
        with ``abutils`` will be used.

    id_key : str, default=None
        Key to retrieve the sequence ID. If not provided or missing, ``Sequence.id`` is used.

    sequence_key : str, default=None
        Key to retrieve the sequence. If not provided or missing, ``Sequence.sequence`` is used.

    as_dict : bool, default=False
        If ``True``, return clustering results as a ``dict`` rather than the TSV output file.
        object. the ``dict`` is of the format:
            {"cluster1_name": {"centroid_id": cluster1_centroid_id,
                               "seq_ids": [seq1_id, seq2_id, seq3_id, ...]},
             "cluster2_name": {"centroid_id": cluster2_centroid_id,
                               "seq_ids": [seq4_id, seq5_id, seq6_id, ...]},
            }

    quiet : bool, default=False
        If ``True``, suppresses all output from the clustering algorithm.

    debug : bool, default=False
        If ``True``, prints standard output and standard error from ``mmseqs``.
        Default is ``False``.


    Returns
    -------
    clusters : Path to the TSV output file from ``mmseqs`` or a ``dict`` of cluster info.

    Raises
    ------
    ValueError
        If neither or both of ``fasta_file`` and ``db_path`` are provided.

    """
    # Validate input: exactly one of fasta_file or db_path must be provided
    if fasta_file is None and db_path is None:
        raise ValueError("Either 'fasta_file' or 'db_path' must be provided.")
    if fasta_file is not None and db_path is not None:
        raise ValueError("Only one of 'fasta_file' or 'db_path' should be provided, not both.")

    # Track whether we created the database (for cleanup)
    db_is_temporary = db_path is None

    # output files
    if db_is_temporary:
        db_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, prefix="DB_").name
    else:
        db_file = db_path
    clu_file = tempfile.NamedTemporaryFile(
        dir=temp_dir, delete=False, prefix="CLU_"
    ).name
    tsv_file = tempfile.NamedTemporaryFile(
        dir=temp_dir, delete=False, prefix="TSV_"
    ).name
    # get the mmseqs binary
    if mmseqs_bin is None:
        mmseqs_bin = get_binary_path("mmseqs")
    # build the mmseqs DB (only if fasta_file is provided)
    if db_is_temporary:
        db_cmd = f"{mmseqs_bin} createdb {fasta_file} {db_file}"
        p = sp.Popen(db_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        stdout, stderr = p.communicate()
        if debug:
            print("STDOUT:", stdout)
            print("")
            print("STDERR:", stderr)
    # do the clustering
    cluster_cmd = f"{mmseqs_bin} cluster"
    cluster_cmd += f" {db_file} {clu_file} {temp_dir}"
    cluster_cmd += f" --min-seq-id {threshold}"
    cluster_cmd += f" --cluster-mode {cluster_mode}"
    cluster_cmd += f" --cov-mode {cov_mode}"
    cluster_cmd += f" -c {coverage}"
    cluster_cmd += f" --alignment-mode {alignment_mode}"
    cluster_cmd += f" --seq-id-mode {seq_id_mode}"
    if threads is not None:
        cluster_cmd += f" --threads {threads}"
    p = sp.Popen(cluster_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        print("STDOUT:", stdout)
        print("")
        print("STDERR:", stderr)
    # generate TSV-formatted output
    tsv_cmd = f"{mmseqs_bin} createtsv"
    tsv_cmd += f" {db_file} {db_file} {clu_file} {tsv_file}"
    p = sp.Popen(tsv_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        print("STDOUT:", stdout)
        print("")
        print("STDERR:", stderr)
    # parse TSV output
    if as_dict:
        cluster_info = {}
        name_dict = {}
        with open(tsv_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                c, s = line.split()
                if c not in name_dict:
                    name = "".join(
                        random.choice(string.ascii_uppercase + string.digits)
                        for _ in range(8)
                    )
                    name_dict[c] = name
                name = name_dict[c]
                if name not in cluster_info:
                    cluster_info[name] = {"centroid_id": c, "seq_ids": []}
                cluster_info[name]["seq_ids"].append(s)
        if not debug:
            os.remove(tsv_file)
            if db_is_temporary:
                os.remove(db_file)
            os.remove(clu_file)
        return cluster_info
    else:
        return tsv_file


def cluster_cdhit(
    fasta_file: str,
    threshold: float = 0.975,
    temp_dir: str = "/tmp",
    cdhit_bin: str = None,
    threads: int = 0,
    as_dict: bool = False,
    quiet: bool = False,
    debug: bool = False,
):
    """
    Clusters sequences using `MMseqs2`_.

    Parameters
    ----------
    fasta_file : string
        Path to a FASTA-formatted file. Required. If you'd like to run ``mmseqs``
        using ``Sequence`` objects as input, use ``cluster(algo="mmseqs")``.

    threshold : float, default=0.975
        Identity threshold for clustering. Must be between 0 and 1.

    temp_dir : str, default="/tmp"
        Path to a directory for temporary storage of clustering files.

    cdhit_bin : str, optional
        Path to a CD-HIT executable. If not provided, the CD-HIT binary bundled
        with ``abutils`` will be used.

    threads : int, default=0
        Number of threads to use for clustering. If not provided, all available CPU cores
        will be used.

    id_key : str, default=None
        Key to retrieve the sequence ID. If not provided or missing, ``Sequence.id`` is used.

    sequence_key : str, default=None
        Key to retrieve the sequence. If not provided or missing, ``Sequence.sequence`` is used.

    as_dict : bool, default=False
        If ``True``, return clustering results as a ``dict`` rather than the TSV output file.
        object. the ``dict`` is of the format:
            {"cluster1_name": {"centroid_id": cluster1_centroid_id,
                               "seq_ids": [seq1_id, seq2_id, seq3_id, ...]},
             "cluster2_name": {"centroid_id": cluster2_centroid_id,
                               "seq_ids": [seq4_id, seq5_id, seq6_id, ...]},
            }

    quiet : bool, default=False
        If ``True``, suppresses all output from the clustering algorithm.

    debug : bool, default=False
        If ``True``, prints standard output and standard error from ``mmseqs``.
        Default is ``False``.


    Returns
    -------
    clusters : Path to the `".clstr"` output file from CD-HIT or a ``dict`` of cluster info.


    .. _CD-HIT
        https://github.com/weizhongli/cdhit/wiki

    """
    # output files
    output_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False).name
    clusters_file = output_file + ".clstr"
    # get the CD-HIT binary
    if cdhit_bin is None:
        cdhit_bin = get_binary_path("cdhit")
        # mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # system = platform.system().lower()
        # machine = platform.machine().lower().replace("x86_64", "amd64")
        # cdhit_bin = os.path.join(mod_dir, f"bin/cdhit_{system}_{machine}")
    # run CD-HIT
    # clustering options are as follows:
    #   - d: length at which to truncate sequence names in the output file,
    #        default 20 (0 to truncate at first space)
    #   - l: length of the throw_away_sequences, default 10
    #   - T: number of threads, default 1 (0 uses all available threads)
    #   - M: max memory (in MB), default 800 (0 uses all available memory)
    wordsize = _get_cdhit_wordsize(threshold)
    cluster_cmd = f"{cdhit_bin} -i {fasta_file} -o {output_file}"
    cluster_cmd += f" -c {threshold} -n {wordsize} -T {threads}"
    cluster_cmd += " -d 0 -l 4 -M 0"
    p = sp.Popen(cluster_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        print("STDOUT:", stdout)
        print("")
        print("STDERR:", stderr)
    if os.stat(clusters_file).st_size == 0:
        if not quiet:
            err = f"WARNING: the CD-HIT output file ({clusters_file}) is empty. "
            err += "Please verify that the input data is valid."
            print(err)
        return {}
    # parse CD-HIT output
    if as_dict:
        cluster_info = {}
        with open(clusters_file) as f:
            raw_clusters = f.read().split("\n>")
        for raw_cluster in raw_clusters:
            if not raw_cluster.strip():
                continue
            cluster_lines = raw_cluster.split("\n")
            centroid_id = None
            cluster_name = cluster_lines[0]
            seq_ids = []
            for cluster_line in cluster_lines[1:]:
                if not cluster_line.strip():
                    continue
                seq_id = cluster_line.split(">")[1].split("...")[0]
                seq_ids.append(seq_id)
                if "*" in cluster_line:
                    centroid_id = seq_id
            if centroid_id is None:
                centroid_id = seq_ids[0]
            cluster_info[cluster_name] = {
                "centroid_id": centroid_id,
                "seq_ids": seq_ids,
            }
        if not debug:
            os.remove(output_file)
            os.remove(clusters_file)
        return cluster_info
    else:
        return clusters_file


def _get_cdhit_wordsize(threshold):
    if threshold > 0.7:
        return 5
    elif threshold > 0.6:
        return 4
    elif threshold > 0.5:
        return 3
    return 2


# -----------------
# MMseqs2 utilities
# -----------------


def create_mmseqs_db(
    sequences: Iterable | str,
    db_path: str,
    id_key: str | None = None,
    seq_key: str | None = None,
    temp_dir: str = "/tmp",
    mmseqs_bin: str | None = None,
    quiet: bool = False,
    debug: bool = False,
) -> str:
    """
    Creates an MMseqs2 database from sequences for reuse across multiple operations.

    Parameters
    ----------
    sequences : iterable or string
        Input sequences in any of the following formats:
            1. list of abutils ``Sequence`` objects
            2. FASTA-formatted string
            3. path to a FASTA-formatted file
            4. list of BioPython ``SeqRecord`` objects
            5. list of lists/tuples, of the format ``[sequence_id, sequence]``
            6. list of strings, with each string being a separate sequence.
        Required.

    db_path : str
        Path where the MMseqs2 database will be created. MMseqs2 creates multiple
        files with this prefix (e.g., db_path, db_path.index, db_path_h, etc.).
        The parent directory must exist. Required.

    id_key : str, optional
        Key to retrieve the sequence ID. If not provided, ``Sequence.id`` is used.

    seq_key : str, optional
        Key to retrieve the sequence. If not provided, ``Sequence.sequence`` is used.

    temp_dir : str, default="/tmp"
        Path to a directory for temporary storage during database creation.

    mmseqs_bin : str, optional
        Path to an MMseqs2 executable. If not provided, the MMseqs2 binary bundled
        with ``abutils`` will be used.

    quiet : bool, default=False
        If ``True``, suppresses all output.

    debug : bool, default=False
        If ``True``, prints standard output and standard error from ``mmseqs``.

    Returns
    -------
    db_path : str
        Path to the created MMseqs2 database (same as input db_path).

    Raises
    ------
    FileNotFoundError
        If the parent directory of db_path does not exist.
    RuntimeError
        If MMseqs2 database creation fails.

    Examples
    --------
    >>> from abutils.tl import create_mmseqs_db, linclust
    >>> db_path = create_mmseqs_db(sequences, "/path/to/my_database")
    >>> # Use the database for multiple linclust runs with different parameters
    >>> linclust(db_path, output_tsv="/path/to/results1.tsv", threshold=0.9)
    >>> linclust(db_path, output_tsv="/path/to/results2.tsv", threshold=0.95)

    """
    # Validate output path
    _validate_output_path(db_path)

    # Get mmseqs binary
    if mmseqs_bin is None:
        mmseqs_bin = get_binary_path("mmseqs")

    # Convert sequences to FASTA
    fasta_file = to_fasta(
        sequences, tempfile_dir=temp_dir, id_key=id_key, sequence_key=seq_key
    )

    # Create MMseqs2 database
    db_cmd = f"{mmseqs_bin} createdb {fasta_file} {db_path}"
    p = sp.Popen(db_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()

    if debug:
        print("STDOUT:", stdout.decode("utf-8"))
        print("")
        print("STDERR:", stderr.decode("utf-8"))

    if p.returncode != 0:
        raise RuntimeError(
            f"MMseqs2 database creation failed with return code {p.returncode}.\n"
            f"Command: {db_cmd}\n"
            f"stderr: {stderr.decode('utf-8')}"
        )

    # Clean up temporary FASTA file
    if not debug and os.path.exists(fasta_file):
        os.remove(fasta_file)

    if not quiet:
        print(f"MMseqs2 database created at: {db_path}")

    return db_path


def linclust(
    sequences: Iterable | str,
    output_tsv: str,
    output_fasta: str | None = None,
    threshold: float = 0.9,
    coverage: float = 0.8,
    cov_mode: str = "0",
    alignment_mode: str = "3",
    seq_id_mode: str = "1",
    kmer_per_seq: int = 0,
    kmer_per_seq_scale: float = 0.0,
    threads: int | None = None,
    temp_dir: str = "/tmp",
    mmseqs_bin: str | None = None,
    id_key: str | None = None,
    seq_key: str | None = None,
    quiet: bool = False,
    debug: bool = False,
) -> dict:
    """
    Performs linear-time clustering using MMseqs2's linclust algorithm.

    ``linclust`` is optimized for large-scale sequence clustering and is significantly
    faster than ``mmseqs cluster`` for very large datasets (millions of sequences).
    Unlike ``cluster_mmseqs()``, this function writes output directly to user-specified
    files rather than returning cluster objects.

    Parameters
    ----------
    sequences : iterable, string, or database path
        Input sequences in any of the following formats:
            1. list of abutils ``Sequence`` objects
            2. FASTA-formatted string
            3. path to a FASTA-formatted file
            4. list of BioPython ``SeqRecord`` objects
            5. list of lists/tuples, of the format ``[sequence_id, sequence]``
            6. list of strings, with each string being a separate sequence
            7. path to a pre-generated MMseqs2 database (created with ``create_mmseqs_db()``)
        Required.

    output_tsv : str
        Path where the TSV output will be written. The TSV file contains two columns:
        cluster representative ID and member ID (one row per cluster member).
        The parent directory must exist. Required.

    output_fasta : str, optional
        Path where representative sequences (cluster centroids) will be written
        in FASTA format. If not provided, representative FASTA is not generated.

    threshold : float, default=0.9
        Minimum sequence identity threshold for clustering. Must be between 0 and 1.

    coverage : float, default=0.8
        Coverage threshold for clustering. Must be between 0 and 1.

    cov_mode : str, default="0"
        Coverage mode. Options are:
            - ``"0"``: bidirectional
            - ``"1"``: target coverage
            - ``"2"``: query coverage
            - ``"3"``: target-in-query length coverage

    alignment_mode : str, default="3"
        Alignment mode. Options are:
            - ``"0"``: automatic
            - ``"1"``: only score and end_pos
            - ``"2"``: also start_pos and cov
            - ``"3"``: also seq.id
            - ``"4"``: only ungapped alignment

    seq_id_mode : str, default="1"
        Sequence ID mode. Options are:
            - ``"0"``: alignment length
            - ``"1"``: shorter sequence length
            - ``"2"``: longer sequence length

    kmer_per_seq : int, default=0
        k-mers per sequence. 0 means use the MMseqs2 default.

    kmer_per_seq_scale : float, default=0.0
        Scale k-mers per sequence based on sequence length. 0.0 disables scaling.

    threads : int, optional
        Number of threads to use. If not provided, MMseqs2 determines automatically.

    temp_dir : str, default="/tmp"
        Path to a directory for temporary storage of intermediate files.

    mmseqs_bin : str, optional
        Path to an MMseqs2 executable. If not provided, the MMseqs2 binary bundled
        with ``abutils`` will be used.

    id_key : str, optional
        Key to retrieve the sequence ID. Only used if sequences is not a database path.

    seq_key : str, optional
        Key to retrieve the sequence. Only used if sequences is not a database path.

    quiet : bool, default=False
        If ``True``, suppresses all output from the clustering algorithm.

    debug : bool, default=False
        If ``True``, prints standard output and standard error from ``mmseqs``,
        and retains intermediate files in temp_dir.

    Returns
    -------
    result : dict
        Dictionary containing:
            - ``"tsv_file"``: Path to the output TSV file
            - ``"fasta_file"``: Path to the representative FASTA file (or None)
            - ``"cluster_count"``: Number of clusters

    Raises
    ------
    FileNotFoundError
        If the parent directory of output files does not exist.
    RuntimeError
        If MMseqs2 linclust fails.
    ValueError
        If sequences input is invalid.

    Examples
    --------
    >>> from abutils.tl import linclust
    >>> # Basic usage with FASTA file
    >>> result = linclust(
    ...     "sequences.fasta",
    ...     output_tsv="clusters.tsv",
    ...     output_fasta="representatives.fasta",
    ...     threshold=0.9
    ... )
    >>> print(f"Created {result['cluster_count']} clusters")
    >>> print(f"TSV output: {result['tsv_file']}")

    >>> # Using pre-generated database for multiple runs
    >>> from abutils.tl import create_mmseqs_db
    >>> db = create_mmseqs_db(sequences, "/data/mydb")
    >>> linclust(db, output_tsv="run1.tsv", threshold=0.9)
    >>> linclust(db, output_tsv="run2.tsv", threshold=0.95)

    See Also
    --------
    cluster_mmseqs : Standard MMseqs2 clustering (returns Clusters objects)
    create_mmseqs_db : Create reusable MMseqs2 database
    cluster : Auto-selecting clustering with multiple backend support

    """
    # Validate output paths
    _validate_output_path(output_tsv)
    if output_fasta:
        _validate_output_path(output_fasta)

    # Get mmseqs binary
    if mmseqs_bin is None:
        mmseqs_bin = get_binary_path("mmseqs")

    # Track temporary files for cleanup
    temp_files = []
    db_is_temporary = False
    fasta_is_temporary = False

    try:
        # Determine input type and prepare database
        if _is_mmseqs_db(sequences):
            # Use pre-generated database directly
            db_path = sequences
        else:
            # Convert sequences to FASTA and create temporary database
            fasta_file = to_fasta(
                sequences, tempfile_dir=temp_dir, id_key=id_key, sequence_key=seq_key
            )
            fasta_is_temporary = True

            db_path = tempfile.NamedTemporaryFile(
                dir=temp_dir, delete=False, prefix="linclust_db_"
            ).name
            db_is_temporary = True

            # Create database
            db_cmd = f"{mmseqs_bin} createdb {fasta_file} {db_path}"
            p = sp.Popen(db_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            stdout, stderr = p.communicate()

            if debug:
                print("CREATEDB STDOUT:", stdout.decode("utf-8"))
                print("CREATEDB STDERR:", stderr.decode("utf-8"))

            if p.returncode != 0:
                raise RuntimeError(
                    f"MMseqs2 createdb failed with return code {p.returncode}.\n"
                    f"stderr: {stderr.decode('utf-8')}"
                )

        # Create temporary cluster output file
        clu_path = tempfile.NamedTemporaryFile(
            dir=temp_dir, delete=False, prefix="linclust_clu_"
        ).name
        temp_files.append(clu_path)

        # Run linclust
        linclust_cmd = f"{mmseqs_bin} linclust"
        linclust_cmd += f" {db_path} {clu_path} {temp_dir}"
        linclust_cmd += f" --min-seq-id {threshold}"
        linclust_cmd += f" -c {coverage}"
        linclust_cmd += f" --cov-mode {cov_mode}"
        linclust_cmd += f" --alignment-mode {alignment_mode}"
        linclust_cmd += f" --seq-id-mode {seq_id_mode}"
        if kmer_per_seq > 0:
            linclust_cmd += f" --kmer-per-seq {kmer_per_seq}"
        if kmer_per_seq_scale > 0:
            linclust_cmd += f" --kmer-per-seq-scale {kmer_per_seq_scale}"
        if threads is not None:
            linclust_cmd += f" --threads {threads}"

        p = sp.Popen(linclust_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        stdout, stderr = p.communicate()

        if debug:
            print("LINCLUST STDOUT:", stdout.decode("utf-8"))
            print("LINCLUST STDERR:", stderr.decode("utf-8"))

        if p.returncode != 0:
            raise RuntimeError(
                f"MMseqs2 linclust failed with return code {p.returncode}.\n"
                f"stderr: {stderr.decode('utf-8')}"
            )

        # Generate TSV output
        tsv_cmd = f"{mmseqs_bin} createtsv"
        tsv_cmd += f" {db_path} {db_path} {clu_path} {output_tsv}"

        p = sp.Popen(tsv_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        stdout, stderr = p.communicate()

        if debug:
            print("CREATETSV STDOUT:", stdout.decode("utf-8"))
            print("CREATETSV STDERR:", stderr.decode("utf-8"))

        if p.returncode != 0:
            raise RuntimeError(
                f"MMseqs2 createtsv failed with return code {p.returncode}.\n"
                f"stderr: {stderr.decode('utf-8')}"
            )

        # Generate representative FASTA if requested
        if output_fasta:
            # Create representative sequence database
            repseq_path = tempfile.NamedTemporaryFile(
                dir=temp_dir, delete=False, prefix="linclust_repseq_"
            ).name
            temp_files.append(repseq_path)

            repseq_cmd = f"{mmseqs_bin} result2repseq {db_path} {clu_path} {repseq_path}"
            p = sp.Popen(repseq_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            stdout, stderr = p.communicate()

            if debug:
                print("RESULT2REPSEQ STDOUT:", stdout.decode("utf-8"))
                print("RESULT2REPSEQ STDERR:", stderr.decode("utf-8"))

            if p.returncode != 0:
                raise RuntimeError(
                    f"MMseqs2 result2repseq failed with return code {p.returncode}.\n"
                    f"stderr: {stderr.decode('utf-8')}"
                )

            # Convert to FASTA
            fasta_cmd = f"{mmseqs_bin} result2flat {db_path} {db_path} {repseq_path} {output_fasta} --use-fasta-header"
            p = sp.Popen(fasta_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
            stdout, stderr = p.communicate()

            if debug:
                print("RESULT2FLAT STDOUT:", stdout.decode("utf-8"))
                print("RESULT2FLAT STDERR:", stderr.decode("utf-8"))

            if p.returncode != 0:
                raise RuntimeError(
                    f"MMseqs2 result2flat failed with return code {p.returncode}.\n"
                    f"stderr: {stderr.decode('utf-8')}"
                )

        # Count clusters from TSV
        cluster_reps = set()
        with open(output_tsv) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    rep_id = line.strip().split("\t")[0]
                    cluster_reps.add(rep_id)
        cluster_count = len(cluster_reps)

        if not quiet:
            print(f"linclust: clustered sequences into {cluster_count} clusters")
            print(f"TSV output: {output_tsv}")
            if output_fasta:
                print(f"Representative FASTA: {output_fasta}")

        return {
            "tsv_file": output_tsv,
            "fasta_file": output_fasta,
            "cluster_count": cluster_count,
        }

    finally:
        # Clean up temporary files
        if not debug:
            # Clean up cluster output files
            for f in temp_files:
                _cleanup_mmseqs_db(f)

            # Clean up temporary database
            if db_is_temporary:
                _cleanup_mmseqs_db(db_path)

            # Clean up temporary FASTA
            if fasta_is_temporary and os.path.exists(fasta_file):
                os.remove(fasta_file)


# def cluster(seqs, threshold=0.975, out_file=None, temp_dir=None, make_db=True, method='cdhit',
#             quiet=False, threads=0, return_just_seq_ids=False, max_memory=800, retries=5, debug=False):
#     '''
#     Perform sequence clustering with CD-HIT.

#     Args:

#         seqs (list or filepath): An iterable of sequences, in any format that ``abutils.utils.sequence.Sequence``
#             can handle. Alternatively, ``seqs`` can be the path to a FASTA-formatted file of sequences. If a path
#             is provided, ``make_db`` is set to ``False`` and the input file will not be deleted following clustering.

#         threshold (float): Clustering identity threshold. Default is ``0.975``.

#         out_file (str): Path to the clustering output file. Default is to use
#             ``tempfile.NamedTemporaryFile`` to generate an output file name.

#         temp_dir (str): Path to the temporary directory. If not provided, ``'/tmp'`` is used.

#         make_db (bool): Whether to build a SQlite database of sequence information. Required
#             if you want to calculate consensus/centroid sequences for the resulting
#             clusters or if you need to access the clustered sequences (not just sequence IDs)
#             Default is ``True``.

#         method (str): Doesn't do anything at this point, stubbed to allow expansion of sequencing tools
#             beyond CD-HIT (UCLUST, etc).

#         quiet (bool): If ``True``, surpresses printing of output/progress info. Default is ``False``.

#         threads (int): Number of threads (CPU cores) to be used for clustering. Default is ``0``,
#             which results in all available cores being used.

#         return_just_seq_ids (bool): If ``True``, will return a 2D list of sequence IDs
#             (a list containing a list of sequence IDs for each cluster) rather than returning a
#             list of ``Cluster`` objects.

#         max_memory (int): Max memory (in MB) for CD-HIT. Will be passed directly to CD-HIT through
#             the ``-M`` runtime option. Default is ``800``.

#         retries (int): Occasionally (and for no discernable reason), the CD-HIT output files are empty.
#             If this occurs, ``cdhit`` will attempt to retry running CD-HIT until the output files are
#             not empty. Default is 5.

#         debug (bool): If ``True``, print standard output and standard error from CD-HIT. Default is ``False``.

#     Returns:

#         list: A ``CDHITResult`` object or, if ``return_just_seq_ids`` is ``True``, a 2D list of sequence IDs.
#     '''
#     if make_db:
#         ofile, cfile, seq_db, db_path = cdhit(seqs, out_file=out_file, temp_dir=temp_dir,
#                                               threshold=threshold, make_db=True, quiet=quiet,
#                                               threads=threads, max_memory=max_memory, debug=debug)
#         return parse_clusters(ofile, cfile, seq_db=seq_db, db_path=db_path, return_just_seq_ids=return_just_seq_ids)
#     else:
#         if isinstance(seqs, (list, tuple)):
#             seqs = [Sequence(s) for s in seqs]
#             seq_dict = {s.id: s for s in seqs}
#         elif os.path.exists(seqs):
#             seq_dict = None
#         ofile, cfile, = cdhit(seqs, out_file=out_file, temp_dir=temp_dir, threads=threads,
#                               threshold=threshold, make_db=False, quiet=quiet,
#                               max_memory=max_memory, debug=debug)
#         return parse_clusters(ofile, cfile, seq_dict=seq_dict, return_just_seq_ids=return_just_seq_ids)


# def cdhit(seqs, out_file=None, temp_dir=None, threshold=0.975, make_db=True,
#           quiet=False, threads=0, max_memory=800, retries=5, debug=False):
#     '''
#     Run CD-HIT.

#     Args:

#         seqs (list or filepath): An iterable of sequences, in any format that `abutils.utils.sequence.Sequence()`
#             can handle. Alternatively, ``seqs`` can be the path to a FASTA-formatted file of sequences. If a file path
#             is provided, ``make_db`` is set to ``False`` and the input file will not be deleted following clustering.

#         threshold (float): Clustering identity threshold. Default is `0.975`.

#         out_file (str): Path to the clustering output file. Default is to use
#             `tempfile.NamedTemporaryFile` to generate an output file name.

#         temp_dir (str): Path to the temporary directory. If not provided, the default
#             temporary directory is used (typically ``/tmp`` on MacOS and Linux).

#         make_db (bool): Whether to build a SQlite database of sequence information. Required
#             if you want to calculate consensus/centroid sequences for the resulting
#             clusters or if you need to access the clustered sequences (not just sequence IDs).
#             Default is `True`.

#         quiet (bool): If `True`, surpresses printing of output/progress info. Default is `False`.

#         threads (int): Number of threads (CPU cores) to be used for clustering. Default is `0`,
#             which results in all available cores being used.

#         max_memory (int): Max memory (in MB) for CD-HIT. Will be passed directly to CD-HIT through
#             the `-M` runtime option. Default is `800`.

#         retries (int): Occasionally (and for no discernable reason), the CD-HIT output files are empty.
#             If this occurs, ``cdhit`` will attempt to retry running CD-HIT until the output files are
#             not empty. Default is 5.

#         debug (bool): If `True`, print standard output and standard error from CD-HIT. Default is `False`.

#     Returns:

#         If `make_db` is `True`, returns the CD-HIT output file path, the CD-HIT cluster file path,
#             a `sqlite3` database connection object, and the database path. If `make_db` is `False`, only the
#             CD-HIT output file path and CD-HIT cluster file path are returned.
#     '''
#     start_time = time.time()
#     delete_input = False if debug else True
#     if not quiet:
#         print('CD-HIT: clustering {} seqeunces'.format(len(seqs)))
#     if out_file is None:
#         out_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
#         out_file.close()
#         ofile = out_file.name
#     else:
#         ofile = os.path.expanduser(out_file)
#     cfile = ofile + '.clstr'
#     with open(ofile, 'w') as f: f.write('')
#     with open(cfile, 'w') as f: f.write('')
#     if isinstance(seqs, (list, tuple)):
#         seqs = [Sequence(s) for s in seqs]
#         ifile = _make_cdhit_input(seqs, temp_dir)
#     elif os.path.exists(seqs):
#         ifile = seqs
#         delete_input = False
#     word_size = _get_cdhit_wordsize(threshold)
#     cdhit_cmd = 'cdhit -i {} -o {} -c {} -n {} -d 0 -T {} -M {}'.format(ifile,
#                                                                         ofile,
#                                                                         threshold,
#                                                                         word_size,
#                                                                         threads,
#                                                                         max_memory)
#     while not all([os.path.getsize(cfile), os.path.getsize(cfile)]):
#         cluster = sp.Popen(cdhit_cmd,
#                         shell=True,
#                         stdout=sp.PIPE,
#                         stderr=sp.PIPE)
#         stdout, stderr = cluster.communicate()
#         if not retries:
#             break
#         retries -= 1
#     end_time = time.time()
#     if debug:
#         print(stdout)
#         print(stderr)
#     if delete_input:
#         os.unlink(ifile)
#     if not quiet:
#         print('CD-HIT: clustering took {:.2f} seconds'.format(end_time - start_time))
#     if make_db:
#         if not quiet:
#             print('CD-HIT: building a SQLite3 database')
#         seq_db, db_path = _build_seq_db(seqs, direc=temp_dir)
#         return ofile, cfile, seq_db, db_path
#     return ofile, cfile


# def parse_clusters(out_file, clust_file, seq_db=None, db_path=None, seq_dict=None, return_just_seq_ids=False):
#     '''
#     Parses CD-HIT output.

#     Args:

#         out_file (str): Path to the CD-HIT output file. Required.

#         clust_file (str): Path to the CD-HIT cluster file. Required.

#         seq_db (sqlite.Connection): SQLite3 `Connection` object. Default is ``None``. If not provided and
#             ``return_just_seq_ids`` is False, the returned ``Cluster`` objects will not contain any sequence
#             information beyond the sequence ID.

#         db_path (str): Path to a SQLite3 database file. Default is ``None``. Must be provided if
#             ``seq_db`` is also provided.

#         seq_dict (dict): A ``dict`` mapping sequence IDs to ``abutils.core.sequence.Sequence`` objects. Default
#             is ``None``. Typically used when a relatively small number of sequences are being clustered and
#             creating a ``sqlite3`` database would be overkill.

#         return_just_seq_ids (bool): If ``True``, will return a 2D list of sequence IDs
#             (a list containing a list of sequence IDs for each cluster) rather than returning a
#             list of ``Cluster`` objects.

#     Returns:

#         A ``CDHITResult`` object or, if ``return_just_seq_ids`` is ``True``, a 2D list of sequence IDs.
#     '''
#     raw_clusters = [c.split('\n') for c in open(clust_file, 'r').read().split('\n>')]
#     if return_just_seq_ids:
#         ids = []
#         for rc in raw_clusters:
#             _ids = []
#             for c in rc[1:]:
#                 if c:
#                     _ids.append(c.split()[2][1:-3])
#             ids.append(_ids)
#         os.unlink(out_file)
#         os.unlink(clust_file)
#         return ids
#     os.unlink(out_file)
#     os.unlink(clust_file)
#     clusters = [Cluster(rc, seq_db, db_path, seq_dict) for rc in raw_clusters]
#     return CDHITResult(clusters, seq_db=seq_db, db_path=db_path, seq_dict=seq_dict)


# def _make_cdhit_input(seqs, temp_dir):
#     ifile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
#     ifile.close()
#     fastas = [s.fasta for s in seqs]
#     with open(ifile.name, 'w') as f:
#         f.write('\n'.join(fastas))
#     return ifile.name


# def _get_cdhit_wordsize(threshold):
#     if threshold > 0.7:
#         return 5
#     elif threshold > 0.6:
#         return 4
#     elif threshold > 0.5:
#         return 3
#     return 2


# def _build_seq_db(seqs, direc=None):
#     # '''
#     # Builds a SQLite3 database of sequences.

#     # Inputs are a list of Sequence objects and an optional directory to store the database.
#     # If ::direc:: is not provided, '/tmp' will be used.

#     # Returns a SQLite3 connection object and the database path.
#     # '''
#     direc = direc if direc is not None else '/tmp'
#     db_name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
#     db_path = os.path.join(direc, db_name)
#     conn = sqlite3.connect(db_path)
#     c = conn.cursor()
#     create_cmd = '''CREATE TABLE seqs (id text, sequence text)'''
#     insert_cmd = 'INSERT INTO seqs VALUES (?,?)'
#     c.execute('DROP TABLE IF EXISTS seqs')
#     c.execute(create_cmd)
#     c.executemany(insert_cmd, [(str(s.id), str(s.sequence)) for s in seqs])
#     c.execute('CREATE INDEX seq_index ON seqs (id)')
#     return c, db_path
