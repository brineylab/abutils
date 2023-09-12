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


__all__ = ["cluster", "cluster_vsearch", "cluster_mmseqs", "cluster_cdhit"]


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
    cdhit_bin: str = None,
    id_key: Optional[str] = None,
    seq_key: Optional[str] = None,
    strand: str = "plus",
    as_dict: bool = False,
    debug: bool = False,
) -> Union[dict, Clusters]:
    """
    Clusters sequences using `CD-HIT`_, `VSEARCH`_ or `MMseqs2`_. By default, sequences will
    be clustered with VSEARCH if there are fewer than 10,000 nucleotide sequences, with
    CD-HIT if there are fewer than 10,000 amino acide sequences, and with MMseqs2 if there
    are more than 10,000 sequences (nucleotide or amino acid). These defaults can be
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
        Algorithm to be used for clustering. Options are ``"vsearch"``, ``"mmseqs"``, ``"cdhit"``,
        or ``"auto"``. By default (``"auto"``), VSEARCH will be used if there are fewer than 10,000
        nucleotide sequences, CD-HIT will be used if there are fewer than 10,000 amino acid sequences,
        and MMseqs2 will be used for 10,000 sequences or more. Providing ``"vsearch"``, ``"cdhit"``,
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

    cdhit_bin : str, optional
        Path to a CD-HIT executable. If not provided, the CD-HIT binary bundled
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
    .. _CD-HIT
       https://github.com/weizhongli/cdhit/wiki

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
            temp_dir=temp_dir,
            mmseqs_bin=mmseqs_bin,
            as_dict=True,
            debug=debug,
        )
    elif algo in ["cdhit", "cd_hit", "cd-hit"]:
        cluster_dict = cluster_cdhit(
            fasta_file=fasta_file,
            threshold=threshold,
            temp_dir=temp_dir,
            cdhit_bin=cdhit_bin,
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
        system = platform.system().lower()
        machine = platform.machine().lower().replace("x86_64", "amd64")
        vsearch_bin = os.path.join(mod_dir, f"bin/vsearch_{system}_{machine}")
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
    tsv_file = tempfile.NamedTemporaryFile(
        dir=temp_dir, delete=False, prefix="TSV_"
    ).name
    # get the mmseqs binary
    if mmseqs_bin is None:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        system = platform.system().lower()
        machine = platform.machine().lower().replace("x86_64", "amd64")
        mmseqs_bin = os.path.join(mod_dir, f"bin/mmseqs_{system}_{machine}")
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


def cluster_cdhit(
    fasta_file: str,
    threshold: float = 0.975,
    temp_dir: str = "/tmp",
    cdhit_bin: str = None,
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

    cdhit_bin : str, optional
        Path to a CD-HIT executable. If not provided, the CD-HIT binary bundled
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
    clusters : Path to the `".clstr"` output file from CD-HIT or a ``dict`` of cluster info.


    .. _CD-HIT
        https://github.com/weizhongli/cdhit/wiki

    """
    # output files
    output_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False).name
    clusters_file = output_file + ".clstr"
    # get the CD-HIT binary
    if cdhit_bin is None:
        mod_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        system = platform.system().lower()
        machine = platform.machine().lower().replace("x86_64", "amd64")
        cdhit_bin = os.path.join(mod_dir, f"bin/cdhit_{system}_{machine}")
    # run CD-HIT
    # clustering options are as follows:
    #   - d: length at which to truncate sequence names in the output file,
    #        default 20 (0 to truncate at first space)
    #   - l: length of the throw_away_sequences, default 10
    #   - T: number of threads, default 1 (0 uses all available threads)
    #   - M: max memory (in MB), default 800 (0 uses all available memory)
    wordsize = _get_cdhit_wordsize(threshold)
    cluster_cmd = f"{cdhit_bin} -i {fasta_file} -o {output_file}"
    cluster_cmd += f" -c {threshold} -n {wordsize}"
    cluster_cmd += " -d 0 -l 4 -T 0 -M 0"
    p = sp.Popen(cluster_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    if debug:
        print("STDOUT:", stdout)
        print("")
        print("STDERR:", stderr)
    if os.stat(clusters_file).st_size == 0:
        err = f"WARNING: the CD-HIT output file ({clusters_file}) is empty. "
        err += "Please verify that the input data is valid."
        print(err)
        return None
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
