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


import os
import platform
import random
import string
import subprocess as sp
import sys
import tempfile
from typing import Iterable, Optional, Union

from Bio.Align import AlignInfo

from ..utils.alignment import mafft
from ..utils.decorators import lazy_property
from ..utils.pipeline import make_dir

from ..core.sequence import Sequence, read_fasta, to_fasta


class Cluster:
    """
    Docstring for Cluster.
    """

    def __init__(self, name, sequences, centroid=None):
        self.name = name
        self.sequences = sequences
        self.centroid = centroid

    def __iter__(self):
        for s in self.sequences:
            yield s

    @property
    def size(self):
        return len(self.sequences)

    @property
    def seq_ids(self):
        return [s.id for s in self.sequences]

    @lazy_property
    def consensus(self):
        return self._make_consensus()

    def _make_consensus(self):
        if len(self.sequences) == 1:
            return self.sequences[0]
        aln = mafft(self.sequences)
        if aln is None:
            print("ERROR: Failed to generate an alignmnet for a consensus sequence.")
            return None
        summary_align = AlignInfo.SummaryInfo(aln)
        consensus = summary_align.gap_consensus(threshold=0.51, ambiguous="n")
        consensus_string = str(consensus).replace("-", "")
        consensus_seq = Sequence(consensus_string.upper())
        return consensus_seq


class Clusters:
    """
    Docstring for Clusters.
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
                seqs = cdata["seqs"]
                centroid = cdata["centroid"]
                clusters.append(Cluster(name, seqs, centroid=centroid))
            return clusters
        elif clusts is None:
            return []
        else:
            return clusts

    def add(self, cluster):
        self._clusters.append(cluster)


def cluster(
    sequences: Union[Iterable, str],
    threshold: float = 0.975,
    algo: str = "auto",
    temp_dir: str = "/tmp",
    iddef: int = 0,
    vsearch_bin: str = None,
    mmseqs_bin: str = None,
    id_key: Optional[str] = None,
    seq_key: Optional[str] = None,
    strand: str = "plus",
    as_dict: bool = False,
    debug: bool = False,
) -> Union[dict, Clusters]:
    """
    Clusters sequences using `VSEARCH`_ or `MMseqs2`_. By default, sequences will
    be clustered with VSEARCH if there are fewer than 10,000 sequences and
    with MMseqs2 if there are more than 10,000 sequences. These defaults can be
    overridden with `algo`.

    Parameters
    ----------
    sequences : iterable or string
        Input sequences in any of the following formats:
            1. list of abutils ``Sequence`` objects
            2. FASTA-formatted string
            3. path to a FASTA-formatted file
            4. list of BioPython ``SeqRecord`` objects
            5. list of lists/tuples, of the format ``[sequence_id, sequence]``
        Required.

    threshold : float, default=0.975
        Identity threshold for clustering. Must be between 0 and 1.

    algo : float, default="auto"
        Algorithm to be used for clustering. Options are ``"vsearch"``, ``"mmseqs"``,
        or ``"auto"``. By default (``"auto"``), VSEARCH will be used if there are fewer than 10,000
        sequences and MMseqs2 will be used for 10,000 sequences or more. Providing ``"vsearch"``
        or ``"mmseqs"`` will force the use of the desired clustering algorithm regardless of the
        number or sequences to be clustered.

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

    mmseqs_bin : str, optional
        Path to a MMseqs2 executable. If not provided, the MMseqs2 binary bundled
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

    debug : bool, default=False
        If ``True``, prints MAFFT's standard output and standard error.
        Default is ``False``.


    Returns
    -------
    clusters : ``Clusters`` or ``dict``


    .. _VSEARCH
        https://github.com/torognes/vsearch
        x
    .. _MMseqs2
        https://github.com/soedinglab/MMseqs2

    """
    # check input data to get the number of sequences
    fasta_file = to_fasta(
        sequences, tempfile_dir=temp_dir, id_key=id_key, sequence_key=seq_key
    )
    seqs = read_fasta(fasta_file)
    seq_dict = {s.id: s for s in seqs}
    # select the clustering algo
    algo = algo.lower()
    if algo == "auto" and len(seqs) < 10000:
        algo = "vsearch"
    elif algo == "auto" and len(seqs) >= 10000:
        algo = "mmseqs"
    if algo in ["mmseqs", "mmseqs2"]:
        cluster_dict = cluster_mmseqs(
            fasta_file=fasta_file,
            threshold=threshold,
            temp_dir=temp_dir,
            mmseqs_bin=mmseqs_bin,
            as_dict=True,
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
            as_dict=True,
            debug=debug,
        )
    else:
        err = f"\nERROR: Invalid algo option: {algo}."
        err += " Valid choices are: 'vsearch', 'mmseqs', or 'auto'.\n"
        print(err)
        sys.exit()
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
    as_dict: bool = False,
    debug: bool = False,
) -> Union[dict, Clusters]:
    """
    Clusters sequences using `VSEARCH`_.

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

    debug : bool, default=False
        If ``True``, prints standard output and standard error from ``vsearch``.
        Default is ``False``.


    Returns
    -------
    clusters : Path to the UC output file from ``vsearch`` or a ``dict`` of cluster info.


    .. _VSEARCH
        https://github.com/torognes/vsearch

    """
    # output files
    centroid_file = tempfile.NamedTemporaryFile(
        dir=temp_dir, delete=False, prefix="centroids_"
    ).name
    uc_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, prefix="uc_").name
    # get the vsearch binary
    if vsearch_bin is None:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        vsearch_bin = os.path.join(
            mod_dir, "bin/vsearch_{}".format(platform.system().lower())
        )
    # do clustering
    vsearch_cmd = f"{vsearch_bin} --cluster_fast {fasta_file}"
    vsearch_cmd += f" --centroids {centroid_file}"
    vsearch_cmd += f" --clusterout_id"
    vsearch_cmd += f" --uc {uc_file}"
    vsearch_cmd += f" --id {threshold}"
    vsearch_cmd += f" --iddef {iddef}"
    vsearch_cmd += f" --sizeout"
    vsearch_cmd += f" --strand {strand}"
    p = sp.Popen(vsearch_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        print("STDOUT:", stdout.decode("utf-8"))
        print("")
        print("STDERR:", stderr.decode("utf-8"))
    # process output
    if os.stat(uc_file).st_size == 0:
        err = f"WARNING: the VSEARCH output file ({uc_file}) is empty. "
        err += "Please verify that the input data is valid."
        print(err)
        return None
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
        return cluster_info
    else:
        return uc_file


def cluster_mmseqs(
    fasta_file: str,
    threshold: float = 0.975,
    temp_dir: str = "/tmp",
    mmseqs_bin: str = None,
    as_dict: bool = False,
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

    debug : bool, default=False
        If ``True``, prints standard output and standard error from ``mmseqs``.
        Default is ``False``.


    Returns
    -------
    clusters : Path to the TSV output file from ``mmseqs`` or a ``dict`` of cluster info.


    .. _MMseqs2
        https://github.com/soedinglab/MMseqs2

    """
    # output files
    db_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, prefix="DB_").name
    clu_file = tempfile.NamedTemporaryFile(
        dir=temp_dir, delete=False, prefix="CLU_"
    ).name
    tsv_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, prefix="TSV_")
    # get the mmseqs binary
    if mmseqs_bin is None:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        mmseqs_bin = os.path.join(
            mod_dir, "bin/mmseqs_{}".format(platform.system().lower())
        )
    # build the mmseqs DB
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
        return cluster_info
    else:
        return tsv_file
