#!/usr/bin/env python
# filename: clonify.py
#
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


import itertools
import os
from collections import Counter
from typing import Iterable, Optional, Union

import abstar
import fastcluster
import polars as pl
from mnemonic import Mnemonic
from rapidfuzz.distance.Levenshtein import distance as levenshtein_distance
from scipy.cluster.hierarchy import fcluster

from ..core.sequence import Sequence
from ..io import from_polars, make_dir, read_airr, read_parquet
from ..tools.cluster import cluster


def clonify(
    sequences: Iterable[Sequence],
    output_path: Optional[str] = None,
    distance_cutoff: float = 0.32,
    shared_mutation_bonus: float = 0.65,
    length_penalty_multiplier: Union[int, float] = 2,
    group_by_v: bool = True,
    group_by_j: bool = True,
    precluster: bool = False,
    preclustering_threshold: float = 0.65,
    vgene_key: str = "v_gene",
    jgene_key: str = "j_gene",
    cdr3_key: str = "cdr3_aa",
    mutations_key: str = "v_mutations",
    preclustering_key: str = "cdr3_nt",
    mutation_delimiter: str = "|",
    lineage_field: str = "lineage",
    lineage_size_field: str = "lineage_size",
    input_fmt: str = "airr",
    output_fmt: str = "airr",
    temp_directory: Optional[str] = None,
    return_assignment_dict: bool = False,
) -> Union[dict, Iterable[Sequence]]:
    """
    Assigns BCR sequences to clonal lineages using the clonify_ [Briney16]_ algorithm.

    .. seealso::
       | Bryan Briney, Khoa Le, Jiang Zhu, and Dennis R Burton
       | Clonify: unseeded antibody lineage assignment from next-generation sequencing data.
       | *Scientific Reports* 2016. https://doi.org/10.1038/srep23901

    Parameters
    ----------
    sequences : Union[str, Iterable[Sequence]]
        Input sequences in any of the following formats:

            1. list of abutils ``Sequence`` objects
            2. path to a FASTA/Q-formatted file (if `input_fmt` == "fasta" or "fastq")
            3. path to an TSV file of AIRR-formatted annotations (if `input_fmt` == "airr")
            4. path to an Parquet file of AIRR-formatted annotations (if `input_fmt` == "parquet")

        If a FASTA/Q file is provided, ``abstar`` will be used to annotate the sequences
        before lineage assignment.

    output_path: Optional[str], default=None
        Output file path. If not provided, outputs will not be written to file
        and instead, ``Sequence`` objects will be returned. Format of the output file
        is determined by `ouptut_fmt`, and can be either ``"airr"`` (which is TSV-formatted)
        or ``"parquet"``.

    distance_cutoff : float, default=0.32
        Distance threshold for lineage clustering.

    shared_mutation_bonus : float, default=0.65
        Bonus applied for each shared V-gene mutation.

    length_penalty_multiplier : int, default=2
        Multiplier for the CDR3 length penalty. Default is ``2``, resulting in CDR3s that
        differ by ``n`` amino acids being penalized ``n * 2``.

    group_by_v : bool, default=True
        If ``True``, sequences are grouped by V-gene use prior to lineage assignment. This
        option is additive with ``group_by_j``. For example, if ``group_by_v == True`` and
        ``group_by_j == True``, sequences will be grouped by both V-gene and J-gene.

    group_by_j : bool, default=True
        If ``True``, sequences are grouped by J-gene use prior to lineage assignment. This
        option is additive with ``group_by_v``. For example, if ``group_by_v == True`` and
        ``group_by_j == True``, sequences will be grouped by both V-gene and J-gene.

    precluster : bool, default=False
        If ``True``, V/J groups are pre-clustered on the `preclustering_field` sequence,
        which can potentially speed up lineage assignment and reduce memory usage.
        If ``False``, each V/J group is processed in its entirety without pre-clustering.

    preclustering_threshold : float, default=0.65
        Identity threshold for pre-clustering the V/J groups prior to lineage assignment.

    vgene_key : str, default='v_gene'
        Annotation field containing the V-gene name.

    jgene_key : str, default='j_gene'
        Annotation field containing the J-gene name.

    cdr3_key : str, default='cdr3_aa'
        Annotation field containing the CDR3 amino acid sequence.

    muts_key : str, default='v_mutations'
        Annotation field containing the V-gene mutations.

    preclustering_key : str, default='cdr3_nt'
        Annotation field on which to pre-cluster sequences.

    muts_delimiter : str, default='|'
        Delimiter used to separate mutations in the `muts_key` annotation field.

    lineage_field : str, default='lineage'
        Name of the lineage assignment field.

    lineage_size_field : str, default='lineage_size'
        Name of the lineage size field.

    input_fmt : str, default='airr'
        Input format of the sequences. Can be ``"fasta"``, ``"fastq"``, ``"airr"`` or ``"parquet"``.

    output_fmt : str, default='airr'
        Output format of the sequences. Can be ``"airr"`` or ``"parquet"``.

    temp_directory : Optional[str], default=None
        Temporary directory for intermediate files. If not provided, a default
        directory will be created automatically.

    return_assignment_dict : bool, default=False
        If ``True``, a dictionary linking sequence IDs to lineage names will be returned.
        If ``False``, the input ``anndata.AnnData`` object will be returned, with lineage
        annotations included.

    Returns
    -------
    output : ``Iterable[Sequence]`` or ``dict``
        By default (``return_assignment_dict == False``), a list of input sequences are
        returned with two additional annotations: the lineage name (added to `lineage_field`)
        and lineage size (added to `lineage_size_field`). If ``return_assignment_dict == True``,
        a ``dict`` mapping sequence IDs to lineage names is returned.


    .. _clonify: https://github.com/briney/clonify
    """
    # set up file paths
    if output_path is not None:
        output_path = os.path.abspath(output_path)
        make_dir(os.path.dirname(output_path))
    if temp_directory is None:
        if output_path is None:
            temp_directory = "/tmp/.clonify_temp"
        else:
            temp_directory = os.path.join(os.path.dirname(output_path), ".clonify_temp")
    make_dir(temp_directory)

    # process input data
    if isinstance(sequences, str):
        input_fmt = input_fmt.lower()
        if input_fmt in ["fasta", "fastq"]:
            sequences = abstar.run(sequences)
        elif input_fmt == "airr":
            sequences = read_airr(sequences)
        elif input_fmt == "parquet":
            sequences = read_parquet(sequences)
        else:
            raise ValueError(f"Invalid input format: {input_fmt}")
    df = pl.DataFrame(s.annotations for s in sequences)

    # filter DataFrame
    fields = ["sequence_id", vgene_key, jgene_key, cdr3_key, mutations_key]
    if precluster:
        fields.append(preclustering_key)
    filtered_df = df.filter(pl.col("locus") == "IGH")  # just heavy chains
    filtered_df = filtered_df.select(fields)

    # group sequences by V/J genes
    if not group_by_v and not group_by_j:
        group_dfs = [(None, filtered_df)]  # matches polars groupby output
    else:
        group_by = []
        if group_by_v:
            group_by.append(vgene_key)
        if group_by_j:
            group_by.append(jgene_key)
        group_dfs = filtered_df.group_by(group_by)

    # assign lineages
    assign_dict = {}
    mnemo = Mnemonic("english")
    for _, group_df in group_dfs:
        seqs = from_polars(group_df, sequence_key=cdr3_key)
        for s in seqs:
            muts = s[mutations_key]
            s["mutations"] = [] if muts is None else muts.split(mutation_delimiter)
        # preclustering
        if precluster:
            clusters = cluster(
                sequences=seqs,
                seq_key=preclustering_key,
                threshold=preclustering_threshold,
            )
            groups = [c.sequences for c in clusters]
        else:
            groups = [seqs]
        # process preclustered groups
        for group in groups:
            if len(group) == 1:
                # singletons
                seq = group[0]
                assign_dict[seq.id] = "_".join(mnemo.generate(strength=128).split()[:6])
            else:
                # calculate distances
                distances = []
                for s1, s2 in itertools.combinations(group, 2):
                    d = pairwise_distance(
                        s1=s1,
                        s2=s2,
                        vgene_field=vgene_key,
                        jgene_field=jgene_key,
                        cdr3_field=cdr3_key,
                        mutations_field="mutations",
                        shared_mutation_bonus=shared_mutation_bonus,
                        length_penalty_multiplier=length_penalty_multiplier,
                    )
                    distances.append(d)
                # cluster
                linkage_matrix = fastcluster.linkage(
                    distances,
                    method="average",
                    preserve_input=False,
                )
                cluster_list = fcluster(
                    linkage_matrix,
                    distance_cutoff,
                    criterion="distance",
                )
                # rename clusters
                cluster_ids = set(cluster_list)
                cluster_names = {
                    c: "_".join(mnemo.generate(strength=128).split()[:6])
                    for c in cluster_ids
                }
                renamed_clusters = [cluster_names[c] for c in cluster_list]
                # assign sequences
                for seq, name in zip(group, renamed_clusters):
                    assign_dict[seq.id] = name

    # return just the assignment dict if requested
    if return_assignment_dict:
        return assign_dict

    # otherwise, add the lineage name and size to the sequence DataFrame
    lineage_size_dict = Counter(assign_dict.values())
    lineages = [assign_dict.get(s, None) for s in df["sequence_id"]]
    lineage_sizes = [lineage_size_dict.get(lineage, 1) for lineage in lineages]
    df = df.with_columns(
        pl.Series(name=lineage_field, values=lineages),
        pl.Series(name=lineage_size_field, values=lineage_sizes),
    )

    # output
    if output_path is not None:
        if output_fmt == "airr":
            pl.write_csv(df, output_path, separator="\t")
        elif output_fmt == "parquet":
            pl.write_parquet(df, output_path)
    else:
        return from_polars(df)


def pairwise_distance(
    s1: Sequence,
    s2: Sequence,
    shared_mutation_bonus: float = 0.65,
    length_penalty_multiplier: Union[int, float] = 2,
    vgene_field: str = "v_gene",
    jgene_field: str = "j_gene",
    cdr3_field: str = "cdr3_aa",
    mutations_field: str = "mutations",
) -> float:
    """
    Computes length and mutation adjusted Levenshtein distance for a pair of sequences.

    Parameters
    ----------
    s1 : Sequence
        The first input sequence

    s2 : Sequence
        The second input sequence

    shared_mutation_bonus : float, optional
        The bonus for each shared mutation, by default 0.65

    length_penalty_multiplier : Union[int, float], optional
        Used to compute the penalty for differences in CDR3 length. The length
        difference is multiplied by `length_penalty_multiplier`, by default 2

    vgene_field : str, optional
        Name of the field in `s1` and `s2` containing the V-gene name,
        by default ``"v_gene"``

    jgene_field : str, optional
        Name of the field in `s1` and `s2` containing the J-gene name,
        by default ``"j_gene"``

    cdr3_field : str, optional
        Name of the field in `s1` and `s2` containing the CDR3 sequence,
        by default ``"cdr3_aa"``

    mutations_field : str, optional
        Name of the field in `s1` and `s2` containing mutation information,
        by default ``"mutations"`` Note that this field should be a list of
        mutation strings

    Returns
    -------
    distance: float
        The adjusted Levenshtein distance between the two sequences
    """
    germline_penalty = 0
    if s1[vgene_field] != s2[vgene_field]:
        germline_penalty += 10
    if s1[jgene_field] != s2[jgene_field]:
        germline_penalty += 5
    dist = levenshtein_distance(s1[cdr3_field], s2[cdr3_field])
    s1_len = len(s1[cdr3_field])
    s2_len = len(s2[cdr3_field])
    length_penalty = abs(s1_len - s2_len) * length_penalty_multiplier
    length = min(s1_len, s2_len)
    mutation_bonus = (
        len(set(s1[mutations_field]) & set(s2[mutations_field])) * shared_mutation_bonus
    )
    score = (dist + length_penalty - mutation_bonus) / length
    return max(score, 0.001)  # distance values can't be negative


# def clonify(
#     sequences: Iterable[Sequence],
#     distance_cutoff: float = 0.32,
#     shared_mutation_bonus: float = 0.65,
#     length_penalty_multiplier: Union[int, float] = 2,
#     group_by_v: bool = True,
#     group_by_j: bool = True,
#     preclustering: bool = False,
#     preclustering_threshold: float = 0.65,
#     vgene_key: str = "v_gene",
#     jgene_key: str = "j_gene",
#     cdr3_key: str = "cdr3_aa",
#     muts_key: str = "v_mutations",
#     preclustering_key: str = "cdr3_nt",
#     muts_delimiter: str = "|",
#     lineage_field: str = "lineage",
#     lineage_size_field: str = "lineage_size",
#     return_assignment_dict: bool = False,
# ) -> Union[dict, Iterable[Sequence]]:
#     """
#     Assigns BCR sequences to clonal lineages using the clonify_ [Briney16]_ algorithm.

#     .. seealso::
#        | Bryan Briney, Khoa Le, Jiang Zhu, and Dennis R Burton
#        | Clonify: unseeded antibody lineage assignment from next-generation sequencing data.
#        | *Scientific Reports* 2016. https://doi.org/10.1038/srep23901

#     Parameters
#     ----------
#     sequences : Iterable[Sequence]
#         Iterable of ``abutils.core.Sequence`` objects containing annotated sequence data.

#     distance_cutoff : float, default=0.32
#         Distance threshold for lineage clustering.

#     shared_mutation_bonus : float, default=0.65
#         Bonus applied for each shared V-gene mutation.

#     length_penalty_multiplier : int, default=2
#         Multiplier for the CDR3 length penalty. Default is ``2``, resulting in CDR3s that
#         differ by ``n`` amino acids being penalized ``n * 2``.

#     group_by_v : bool, default=True
#         If ``True``, sequences are grouped by V-gene use prior to lineage assignment. This
#         option is additive with ``group_by_j``. For example, if ``group_by_v == True`` and
#         ``group_by_j == True``, sequences will be grouped by both V-gene and J-gene.

#     group_by_j : bool, default=True
#         If ``True``, sequences are grouped by J-gene use prior to lineage assignment. This
#         option is additive with ``group_by_v``. For example, if ``group_by_v == True`` and
#         ``group_by_j == True``, sequences will be grouped by both V-gene and J-gene.

#     preclustering : bool, default=False
#         If ``True``, V/J groups are pre-clustered on the `preclustering_field` sequence,
#         which can potentially speed up lineage assignment and reduce memory usage.
#         If ``False``, each V/J group is processed in its entirety without pre-clustering.

#     preclustering_threshold : float, default=0.65
#         Identity threshold for pre-clustering the V/J groups prior to lineage assignment.

#     vgene_key : str, default='v_gene'
#         Annotation field containing the V-gene name.

#     jgene_key : str, default='j_gene'
#         Annotation field containing the J-gene name.

#     cdr3_key : str, default='cdr3_aa'
#         Annotation field containing the CDR3 amino acid sequence.

#     muts_key : str, default='v_mutations'
#         Annotation field containing the V-gene mutations.

#     preclustering_key : str, default='cdr3_nt'
#         Annotation field on which to pre-cluster sequences.

#     muts_delimiter : str, default='|'
#         Delimiter used to separate mutations in the `muts_key` annotation field.

#     lineage_field : str, default='lineage'
#         Name of the lineage assignment field.

#     lineage_size_field : str, default='lineage_size'
#         Name of the lineage size field.

#     return_assignment_dict : bool, default=False
#         If ``True``, a dictionary linking sequence IDs to lineage names will be returned.
#         If ``False``, the input ``anndata.AnnData`` object will be returned, with lineage
#         annotations included.

#     Returns
#     -------
#     output : ``Iterable[Sequence]`` or ``dict``
#         By default (``return_assignment_dict == False``), a list of input sequences are
#         returned with two additional annotations: the lineage name (added to `lineage_field`)
#         and lineage size (added to `lineage_size_field`). If ``return_assignment_dict == True``,
#         a ``dict`` mapping sequence IDs to lineage names is returned.


#     .. _clonify: https://github.com/briney/clonify
#     """
#     # select the appropriate data fields
#     # if annotation_format.lower() == "airr":
#     #     vgene_key = "v_gene"
#     #     jgene_key = "j_gene"
#     #     cdr3_key = "cdr3_aa"
#     #     muts_key = "v_mutations"
#     # elif annotation_format.lower() == "json":
#     #     vgene_key = "v_gene.gene"
#     #     jgene_key = "j_gene.gene"
#     #     cdr3_key = "cdr3_aa"
#     #     muts_key = "var_muts_nt.muts"
#     # else:
#     #     error = "ERROR: "
#     #     error += 'annotation_format must be either "airr" or "json", '
#     #     err += f"but you provided {annotation_format}"
#     #     print("\n")
#     #     print(error)
#     #     print("\n")
#     #     sys.exit()

#     # setup a list of required fields
#     required_fields = ["v_gene", "j_gene", "cdr3", "mutations"]
#     if preclustering:
#         required_fields.append("preclustering")
#     # group sequences by V/J genes
#     vj_group_dict = {}
#     heavies = [s for s in sequences if s["locus"] == "IGH"]
#     for h in heavies:
#         # build new Sequence objects using just the data we need
#         s = Sequence(h.sequence, id=h.id)
#         s["v_gene"] = nested_dict_lookup(h, vgene_key.split("."))
#         s["j_gene"] = nested_dict_lookup(h, jgene_key.split("."))
#         s["cdr3"] = nested_dict_lookup(h, cdr3_key.split("."))
#         muts = h[muts_key]
#         if any([pd.isnull(muts), muts is None]):
#             s["mutations"] = []
#         elif isinstance(muts, str):
#             muts = muts.split(muts_delimiter)
#         else:
#             muts = list(muts)
#             s["mutations"] = [m for m in muts if m.strip()]
#         if preclustering:
#             s["preclustering"] = nested_dict_lookup(h, preclustering_key.split("."))
#         # skip sequences that don't have all required fields
#         if any([s[v] is None for v in required_fields]):
#             continue
#         # group sequences by VJ gene use
#         vj = ""
#         if group_by_v:
#             vj += s["v_gene"]
#         if group_by_j:
#             vj += s["j_gene"]
#         if vj not in vj_group_dict:
#             vj_group_dict[vj] = []
#         vj_group_dict[vj].append(s)
#     # assign lineages
#     assignment_dict = {}
#     mnemo = Mnemonic("english")
#     for vj_group in vj_group_dict.values():
#         # preclustering
#         if preclustering:
#             seq_dict = {s.id: s for s in vj_group}
#             cluster_seqs = [Sequence(s[preclustering_key], id=s.id) for s in vj_group]
#             clusters = cluster(cluster_seqs, threshold=preclustering_threshold)
#             groups = [[seq_dict[i] for i in c.seq_ids] for c in clusters]
#         else:
#             groups = [
#                 vj_group,
#             ]
#         for group in groups:
#             # singletons
#             if len(group) == 1:
#                 seq = group[0]
#                 assignment_dict[seq.id] = "_".join(
#                     mnemo.generate(strength=128).split()[:6]
#                 )
#             # if not a singleton, cluster
#             else:
#                 # calculate distance matrix
#                 dist_matrix = []
#                 for s1, s2 in itertools.combinations(group, 2):
#                     d = pairwise_distance(
#                         s1, s2, shared_mutation_bonus, length_penalty_multiplier
#                     )
#                     dist_matrix.append(d)
#                 # cluster
#                 linkage_matrix = fc.linkage(
#                     dist_matrix, method="average", preserve_input=False
#                 )
#                 cluster_list = fcluster(
#                     linkage_matrix, distance_cutoff, criterion="distance"
#                 )
#                 # rename clusters
#                 cluster_ids = list(set(cluster_list))
#                 cluster_names = {
#                     c: "_".join(mnemo.generate(strength=128).split()[:6])
#                     for c in cluster_ids
#                 }
#                 renamed_clusters = [cluster_names[c] for c in cluster_list]
#                 # assign sequences
#                 for seq, name in zip(vj_group, renamed_clusters):
#                     assignment_dict[seq.id] = name
#     # return just the assignment dict if requested
#     if return_assignment_dict:
#         return assignment_dict
#     # otherwise, add the lineage name and size to the sequences
#     lineage_size_dict = Counter(assignment_dict.values())
#     for s in sequences:
#         if l := assignment_dict.get(s.id, None) is not None:
#             s[lineage_field] = l
#             s[lineage_size_field] = lineage_size_dict.get(l, 1)
#     return sequences


# def pairwise_distance(
#     s1: Sequence,
#     s2: Sequence,
#     shared_mutation_bonus: float = 0.65,
#     length_penalty_multiplier: Union[int, float] = 2,
#     cdr3_field: str = "cdr3",
#     mutations_field: str = "mutations",
# ) -> float:
#     """
#     Computes length and mutation adjusted Levenshtein distance for a pair of sequences.

#     Parameters
#     ----------
#     s1 : Sequence
#         The first input sequence

#     s2 : Sequence
#         The second input sequence

#     shared_mutation_bonus : float, optional
#         The bonus for each shared mutation, by default 0.65

#     length_penalty_multiplier : Union[int, float], optional
#         Used to compute the penalty for differences in CDR3 length. The length
#         difference is multiplied by `length_penalty_multiplier`, by default 2

#     cdr3_field : str, optional
#         Name of the field in `s1` and `s2` containing the CDR3 sequence,
#         by default ``"cdr3"``

#     mutations_field : str, optional
#         Name of the field in `s1` and `s2` containing mutation information,
#         by default ``"mutations"``

#     Returns
#     -------
#     distance: float
#         The adjusted Levenshtein distance between the two sequences
#     """
#     if len(s1[cdr3_field]) == len(s2[cdr3_field]):
#         dist = sum([i != j for i, j in zip(s1[cdr3_field], s2[cdr3_field])])
#     else:
#         dist = distance(s1[cdr3_field], s2[cdr3_field])
#     length_penalty = (
#         abs(len(s1[cdr3_field]) - len(s2[cdr3_field])) * length_penalty_multiplier
#     )
#     length = min(len(s1[cdr3_field]), len(s2[cdr3_field]))
#     shared_mutations = list(set(s1[mutations_field]) & set(s2[mutations_field]))
#     mutation_bonus = len(shared_mutations) * shared_mutation_bonus
#     score = (dist + length_penalty - mutation_bonus) / length
#     return max(score, 0.001)  # distance values can't be negative