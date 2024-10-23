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
import multiprocessing as mp
import os
import random
import string
from collections import Counter
from typing import Iterable, Optional, Union

import fastcluster
import pandas as pd
import polars as pl
from mnemonic import Mnemonic
from natsort import natsorted
from rapidfuzz.distance.Levenshtein import distance as levenshtein_distance
from scipy.cluster.hierarchy import fcluster
from tqdm.auto import tqdm

from ..core.pair import Pair
from ..core.sequence import Sequence
from ..io import from_polars, make_dir, read_airr, read_csv, read_parquet, to_polars
from ..tools.cluster import cluster
from ..utils.utilities import generate_batches


def clonify(
    sequences: Union[str, Iterable[Sequence]],
    output_path: Optional[str] = None,
    distance_cutoff: float = 0.35,
    shared_mutation_bonus: float = 0.35,
    length_penalty_multiplier: Union[int, float] = 2,
    group_by_v: bool = True,
    group_by_j: bool = True,
    group_by_light_chain_vj: bool = True,
    precluster: bool = False,
    preclustering_threshold: float = 0.65,
    id_key: str = "sequence_id",
    vgene_key: str = "v_gene",
    jgene_key: str = "j_gene",
    cdr3_key: str = "cdr3",
    mutations_key: str = "v_mutations",
    preclustering_key: str = "cdr3",
    mutation_delimiter: str = "|",
    ignore_likely_allelic_variants: bool = False,
    allelic_variant_threshold: float = 0.35,
    min_seqs_for_allelic_variants: int = 200,
    lineage_field: str = "lineage",
    lineage_size_field: str = "lineage_size",
    mnemonic_names: bool = True,
    input_fmt: str = "airr",
    output_fmt: str = "airr",
    temp_directory: Optional[str] = None,
    return_assignment_dict: bool = False,
    batch_size: int = 100000,
    n_processes: int = None,
    verbose: bool = True,
    concise_logging: bool = False,
) -> Union[dict, pl.DataFrame, pd.DataFrame, Iterable[Sequence]]:
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

            1. list of abutils ``Sequence`` or ``Pair`` objects
            2. path to a FASTA/Q-formatted file (if `input_fmt` == "fasta" or "fastq")
            3. path to an TSV file of AIRR-formatted annotations (if `input_fmt` == "airr")
            4. path to an Parquet file of AIRR-formatted annotations (if `input_fmt` == "parquet")
            5. path to a CSV file of AIRR-formatted annotations (if `input_fmt` == "csv")

        If a FASTA/Q file is provided, ``abstar`` will be used to annotate the sequences
        before lineage assignment.

        .. note::
            If a CSV or Parquet file is provided, it can contain annotated data for ``Sequence``
            or ``Pair`` objects. The output file will contain the same fields as the input file,
            with the addition of the lineage annotations.

    output_path: Optional[str], default=None
        Output file path. If not provided, outputs will not be written to file
        and instead, ``Sequence`` or ``Pair`` objects will be returned. Format of the output file
        is determined by `ouptut_fmt`, and can be either ``"airr"`` (which is TSV-formatted),
        ``"csv"``, or ``"parquet"``.

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

    group_by_light_chain_vj : bool, default=True
        If ``True``, sequences are grouped by light chain V-gene and J-gene use prior to
        lineage assignment. This option is additive with ``group_by_v`` and ``group_by_j``.
        For example, if ``group_by_v == True`` and ``group_by_j == True`` and
        ``group_by_light_chain_vj == True``, sequences will be grouped by V-gene, J-gene,
        and light chain V-gene and J-gene. If unpaired sequences are provided, this option
        is ignored.

    .. note::
        To ensure that ``clonify`` behaves in a way that is consistent with user expectations,
        if both ``group_by_v`` and ``group_by_j`` are manually set to ``False``, light chain
        grouping will also be ignored.

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

    cdr3_key : str, default='cdr3'
        Annotation field containing the CDR3 sequence.

    mutations_key : str, default='v_mutations'
        Annotation field containing the V-gene mutations.

    preclustering_key : str, default='cdr3'
        Annotation field on which to pre-cluster sequences.

    mutation_delimiter : str, default='|'
        Delimiter used to separate mutations in the `mutations_key` annotation field.

    ignore_likely_allelic_variants : bool, default=False
        If ``True``, likely allelic variants are ignored when computing shared mutations.

    allelic_variant_threshold : float, default=0.35
        Fraction of sequences in a V-gene that must contain a mutation to be considered an allelic variant.

    min_seqs_for_allelic_variants : int, default=200
        Minimum number of sequences in a V-gene required to consider mutations as allelic variants.

    lineage_field : str, default='lineage'
        Name of the lineage assignment field.

    lineage_size_field : str, default='lineage_size'
        Name of the lineage size field.

    mnemonic_names : bool, default=True
        If ``True``, mnemonic names are used for lineages. If ``False``, lineages are
        named using a random string of alphanumeric characters.

    input_fmt : str, default='airr'
        Input format of the sequences. Can be ``"fasta"``, ``"fastq"``, ``"airr"``,
        ``"csv"``, or ``"parquet"``.

    output_fmt : str, default='airr'
        Output format of the sequences. If `output_path` is provided, the available
        formats are ``"airr"``, ``"csv"``, or ``"parquet"``. If `output_path` is not provided,
        the output format can be ``"polars"``, which returns a ``polars.DataFrame``, or
        ``"pandas"``, which returns a ``pandas.DataFrame``.

    temp_directory : Optional[str], default=None
        Temporary directory for intermediate files. If not provided, a default
        directory will be created automatically.

    batch_size : int, default=100000
        Batch size (per process) for pairwise distance calculations using multiprocessing.

    n_processes : int, default=None
        Number of processes to use for parallelization. If ``None``, the number of processes
        will be set to the number of CPU cores available.

    verbose : bool, default=True
        If ``True``, progress bars and other status updates will be printed to the console.

    concise_logging : bool, default=False
        If ``True``, only minimal progress information will be printed, which can be useful for
        interactive use or when pipelineing. NOT IMPLEMENTED YET.

    Returns
    -------
    output : ``Iterable[Union[Sequence, Pair]]`` or ``polars.DataFrame`` or ``pandas.DataFrame``
        If `output_fmt` is ``"polars"``, a ``polars.DataFrame`` is returned. If `output_fmt` is
        ``"pandas"``, a ``pandas.DataFrame`` is returned. Otherwise, the output is a list of input
        ``Sequence`` or ``Pair`` objects with lineage annotations added. If `output_path` is provided
        and `output_fmt` is one of ``"airr"``, ``"csv"``, or ``"parquet"``, the output is written to
        a file in the specified format, and ``None`` is returned.


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
            import abstar  # import here to avoid circular import

            sequences = abstar.run(sequences)
        elif input_fmt == "airr":
            sequences = read_airr(sequences)
        elif input_fmt == "parquet":
            sequences = read_parquet(sequences)
        elif input_fmt == "csv":
            sequences = read_csv(sequences)
        else:
            raise ValueError(f"Invalid input format: {input_fmt}")
    df = to_polars(sequences)
    is_paired = isinstance(sequences[0], Pair)

    # filter DataFrame
    fields = [id_key, vgene_key, jgene_key, cdr3_key, mutations_key]
    if precluster:
        fields.append(preclustering_key)
    # handle Pair objects
    if is_paired:
        id_key = f"{id_key}:0"
        filtered_df = df.filter(pl.col("locus:0") == "IGH")
        _fields = [f"{f}:0" for f in fields]
        if group_by_light_chain_vj:
            _fields.append(f"{vgene_key}:1")
            _fields.append(f"{jgene_key}:1")
        filtered_df = filtered_df.select(_fields).rename(lambda c: c.replace(":0", ""))
    else:
        filtered_df = df.filter(pl.col("locus") == "IGH")
        filtered_df = filtered_df.select(fields)

    # split mutations string into a list
    mutation_lists = [
        [] if m == "" else m.split(mutation_delimiter)
        for m in filtered_df[mutations_key]
    ]
    filtered_df = filtered_df.with_columns(
        pl.Series(name="mutations", values=mutation_lists),
    )

    # identify mutations associated with potential allelic variants
    unique_vgenes = filtered_df[vgene_key].unique().to_list()
    likely_allelic_variants = {}
    if ignore_likely_allelic_variants:
        # for each V-gene, find mutations that are present
        # at a frequency above allelic_variant_threshold
        if verbose:
            print("- identifying mutations that are likely allelic variants...")
        for v in unique_vgenes:
            v_muts = []
            v_df = filtered_df.filter(pl.col(vgene_key) == v)
            if v_df.shape[0] < min_seqs_for_allelic_variants:
                likely_allelic_variants[v] = []
                continue
            allele_threshold = v_df.shape[0] * allelic_variant_threshold
            for muts in v_df["mutations"]:
                v_muts.extend(muts)
            mut_counts = Counter(v_muts)
            likely_allelic_variants[v] = [
                m for m, c in mut_counts.items() if c >= allele_threshold
            ]
        if verbose:
            for v, muts in likely_allelic_variants.items():
                if muts:
                    print(f"    {v}: {', '.join(natsorted(muts))}")

    else:
        # just make an empty list of "allelic mutations" for each V-gene
        for v in unique_vgenes:
            likely_allelic_variants[v] = []

    # group sequences by V/J genes
    if not group_by_v and not group_by_j:
        group_dfs = [filtered_df]
    else:
        group_by = []
        group_by_list = []
        if group_by_v:
            group_by.append(vgene_key)
            group_by_list.append("V gene")
        if group_by_j:
            group_by.append(jgene_key)
            group_by_list.append("J gene")
        if is_paired and group_by_light_chain_vj:
            group_by.append(f"{vgene_key}:1")
            group_by.append(f"{jgene_key}:1")
            group_by_list.append("Light chain V/J genes")
        if verbose:
            print(f"- grouping by {' and '.join(group_by_list)}")
        grouped = filtered_df.group_by(group_by)
        group_dfs = [g[1] for g in grouped]
        # group_dfs = sorted(group_dfs, key=lambda x: x.shape[0], reverse=True)

    # preclustering
    sequence_groups = []
    if precluster:
        if verbose:
            print("- preclustering")
        for group_df in group_dfs:
            seqs = from_polars(group_df, sequence_key=cdr3_key)
            clusters = cluster(
                sequences=seqs,
                seq_key=preclustering_key,
                threshold=preclustering_threshold,
            )
            sequence_groups.extend([c.sequences for c in clusters])
    else:
        for group_df in group_dfs:
            sequence_groups.append(from_polars(group_df, sequence_key=cdr3_key))

    # configure multiprocessing
    n_processes = n_processes or mp.cpu_count()
    use_mp = True
    if n_processes == 1:
        use_mp = False
    if all([len(sg) < 1000 for sg in sequence_groups]):
        use_mp = False
    if use_mp:
        pool = mp.Pool(processes=n_processes)

    # assign lineages
    mnemo = Mnemonic("english")
    assign_dict = {}
    assign_kwargs = {
        "shared_mutation_bonus": shared_mutation_bonus,
        "length_penalty_multiplier": length_penalty_multiplier,
        "vgene_field": vgene_key,
        "jgene_field": jgene_key,
        "cdr3_field": cdr3_key,
        "mutations_field": "mutations",
    }

    # progress bar
    if verbose:
        print("- assigning lineages:")
        sequence_groups = tqdm(
            sequence_groups,
            position=4,
            # bar_format="{l_bar}{bar:10}| {n_fmt}/{total_fmt} [{elapsed}]"
        )

    # clonify
    for seqs in sequence_groups:
        # no need to process if there's only one sequence
        if len(seqs) == 1:
            if mnemonic_names:
                assign_dict[seqs[0].id] = "_".join(
                    mnemo.generate(strength=128).split()[:8]
                )
            else:
                assign_dict[seqs[0].id] = "".join(
                    random.choices(string.ascii_letters + string.digits, k=16)
                )
            continue

        # distances
        distances = []
        if use_mp:
            # multiprocessing
            async_results = []
            index_iter = itertools.combinations(list(range(len(seqs))), 2)
            for index_batch in generate_batches(index_iter, batch_size):
                async_results.append(
                    pool.apply_async(
                        batch_pairwise_distance,
                        args=(seqs, index_batch),
                        kwds=assign_kwargs,
                    )
                )
            for ar in async_results:
                distances.extend(ar.get())
        else:
            # single processes
            for s1, s2 in itertools.combinations(seqs, 2):
                distances.append(pairwise_distance(s1, s2, **assign_kwargs))

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
        if mnemonic_names:
            cluster_names = {
                c: "_".join(mnemo.generate(strength=128).split()[:8])
                for c in cluster_ids
            }
        else:
            cluster_names = {
                c: "".join(random.choices(string.ascii_letters + string.digits, k=16))
                for c in cluster_ids
            }
        renamed_clusters = [cluster_names[c] for c in cluster_list]
        for seq, name in zip(seqs, renamed_clusters):
            assign_dict[seq.id] = name

    # cleanup multiprocessing (if necessary)
    if use_mp:
        pool.close()
        pool.join()

    # add the lineage name and size to the sequence DataFrame
    lineage_size_dict = Counter(assign_dict.values())
    lineages = [assign_dict.get(s, None) for s in df[id_key]]
    lineage_sizes = [lineage_size_dict.get(lineage, None) for lineage in lineages]
    df = df.with_columns(
        pl.Series(name=lineage_field, values=lineages),
        pl.Series(name=lineage_size_field, values=lineage_sizes),
    )

    # output
    if return_assignment_dict:
        return assign_dict
    if output_path is not None:
        if output_fmt.lower() == "airr":
            pl.write_csv(df, output_path, separator="\t")
        elif output_fmt.lower() == "parquet":
            pl.write_parquet(df, output_path)
        elif output_fmt.lower() == "csv":
            pl.write_csv(df, output_path)
        else:
            raise ValueError
    if output_fmt.lower() == "polars":
        return df
    if output_fmt.lower == "pandas":
        return df.to_pandas()
    return from_polars(df)


def batch_pairwise_distance(sequences, batches, **kwargs):
    distances = []
    for i1, i2 in batches:
        d = pairwise_distance(sequences[i1], sequences[i2], **kwargs)
        distances.append(d)
    return distances


def pairwise_distance(
    s1: Sequence,
    s2: Sequence,
    shared_mutation_bonus: float = 0.35,
    length_penalty_multiplier: Union[int, float] = 2,
    vgene_field: str = "v_gene",
    jgene_field: str = "j_gene",
    cdr3_field: str = "cdr3",
    mutations_field: str = "mutations",
    likely_allelic_variants: Optional[Iterable] = None,
    debug: bool = False,
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

    vgene_field : str, default="v_gene"
        Name of the field in `s1` and `s2` containing the V-gene name,

    jgene_field : str, default="j_gene"
        Name of the field in `s1` and `s2` containing the J-gene name

    cdr3_field : str, default="cdr3"
        Name of the field in `s1` and `s2` containing the CDR3 sequence

    mutations_field : str, default="mutations"
        Name of the field in `s1` and `s2` containing mutation information

    likely_allelic_variants : Optional[Iterable], default=None
        List of mutations that are likely allelic variants. If provided, shared
        mutations that are also in `likely_allelic_variants` will not be counted
        towards the mutation bonus.

    debug : bool, default=False
        If ``True``, debug information will be printed.

    Returns
    -------
    distance: float
        The adjusted Levenshtein distance between the two sequences
    """
    # germline
    germline_penalty = 0
    if s1[vgene_field] != s2[vgene_field]:
        germline_penalty += 10
    if s1[jgene_field] != s2[jgene_field]:
        germline_penalty += 5

    # CDR3 length
    s1_len = len(s1[cdr3_field])
    s2_len = len(s2[cdr3_field])
    length_penalty = abs(s1_len - s2_len) * length_penalty_multiplier
    length = min(s1_len, s2_len)

    # Levenshtein distance
    if s1_len == s2_len:
        # don't allow insertions/deletions if the CDR3s are the same length
        dist = sum([a != b for a, b in zip(s1[cdr3_field], s2[cdr3_field])])
    else:
        dist = levenshtein_distance(s1[cdr3_field], s2[cdr3_field])

    # mutations
    likely_allelic_variants = likely_allelic_variants or []
    mutation_bonus = (
        len(
            set(s1[mutations_field])
            & set(s2[mutations_field]) - set(likely_allelic_variants)
        )
        * shared_mutation_bonus
    )
    score = germline_penalty + ((dist + length_penalty - mutation_bonus) / length)
    return max(score, 0.001)  # distance values can't be negative


# import itertools
# import os
# from collections import Counter
# from typing import Iterable, Optional, Union

# import fastcluster
# import polars as pl
# from mnemonic import Mnemonic
# from rapidfuzz.distance.Levenshtein import distance as levenshtein_distance
# from scipy.cluster.hierarchy import fcluster

# from ..core.sequence import Sequence
# from ..io import from_polars, make_dir, read_airr, read_parquet
# from ..tools.cluster import cluster


# def clonify(
#     sequences: Union[str, Iterable[Sequence]],
#     output_path: Optional[str] = None,
#     distance_cutoff: float = 0.32,
#     shared_mutation_bonus: float = 0.65,
#     length_penalty_multiplier: Union[int, float] = 2,
#     group_by_v: bool = True,
#     group_by_j: bool = True,
#     precluster: bool = False,
#     preclustering_threshold: float = 0.65,
#     vgene_key: str = "v_gene",
#     jgene_key: str = "j_gene",
#     cdr3_key: str = "cdr3_aa",
#     mutations_key: str = "v_mutations",
#     preclustering_key: str = "cdr3_nt",
#     mutation_delimiter: str = "|",
#     lineage_field: str = "lineage",
#     lineage_size_field: str = "lineage_size",
#     input_fmt: str = "airr",
#     output_fmt: str = "airr",
#     temp_directory: Optional[str] = None,
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
#     sequences : Union[str, Iterable[Sequence]]
#         Input sequences in any of the following formats:

#             1. list of abutils ``Sequence`` objects
#             2. path to a FASTA/Q-formatted file (if `input_fmt` == "fasta" or "fastq")
#             3. path to an TSV file of AIRR-formatted annotations (if `input_fmt` == "airr")
#             4. path to an Parquet file of AIRR-formatted annotations (if `input_fmt` == "parquet")

#         If a FASTA/Q file is provided, ``abstar`` will be used to annotate the sequences
#         before lineage assignment.

#     output_path: Optional[str], default=None
#         Output file path. If not provided, outputs will not be written to file
#         and instead, ``Sequence`` objects will be returned. Format of the output file
#         is determined by `ouptut_fmt`, and can be either ``"airr"`` (which is TSV-formatted)
#         or ``"parquet"``.

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

#     precluster : bool, default=False
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

#     input_fmt : str, default='airr'
#         Input format of the sequences. Can be ``"fasta"``, ``"fastq"``, ``"airr"`` or ``"parquet"``.

#     output_fmt : str, default='airr'
#         Output format of the sequences. Can be ``"airr"`` or ``"parquet"``.

#     temp_directory : Optional[str], default=None
#         Temporary directory for intermediate files. If not provided, a default
#         directory will be created automatically.

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
#     # set up file paths
#     if output_path is not None:
#         output_path = os.path.abspath(output_path)
#         make_dir(os.path.dirname(output_path))
#     if temp_directory is None:
#         if output_path is None:
#             temp_directory = "/tmp/.clonify_temp"
#         else:
#             temp_directory = os.path.join(os.path.dirname(output_path), ".clonify_temp")
#     make_dir(temp_directory)

#     # process input data
#     if isinstance(sequences, str):
#         input_fmt = input_fmt.lower()
#         if input_fmt in ["fasta", "fastq"]:
#             import abstar  # here to avoid circular import

#             sequences = abstar.run(sequences)
#         elif input_fmt == "airr":
#             sequences = read_airr(sequences)
#         elif input_fmt == "parquet":
#             sequences = read_parquet(sequences)
#         else:
#             raise ValueError(f"Invalid input format: {input_fmt}")
#     df = pl.DataFrame((s.annotations for s in sequences), infer_schema_length=None)

#     # filter DataFrame
#     fields = ["sequence_id", vgene_key, jgene_key, cdr3_key, mutations_key]
#     if precluster:
#         fields.append(preclustering_key)
#     filtered_df = df.filter(pl.col("locus") == "IGH")  # just heavy chains
#     filtered_df = filtered_df.select(fields)

#     # group sequences by V/J genes
#     if not group_by_v and not group_by_j:
#         group_dfs = [(None, filtered_df)]  # matches polars groupby output
#     else:
#         group_by = []
#         if group_by_v:
#             group_by.append(vgene_key)
#         if group_by_j:
#             group_by.append(jgene_key)
#         group_dfs = filtered_df.group_by(group_by)

#     # assign lineages
#     assign_dict = {}
#     mnemo = Mnemonic("english")
#     for _, group_df in group_dfs:
#         seqs = from_polars(group_df, sequence_key=cdr3_key)
#         for s in seqs:
#             muts = s[mutations_key]
#             s["mutations"] = [] if muts is None else muts.split(mutation_delimiter)
#         # preclustering
#         if precluster:
#             clusters = cluster(
#                 sequences=seqs,
#                 seq_key=preclustering_key,
#                 threshold=preclustering_threshold,
#             )
#             groups = [c.sequences for c in clusters]
#         else:
#             groups = [seqs]
#         # process preclustered groups
#         for group in groups:
#             if len(group) == 1:
#                 # singletons
#                 seq = group[0]
#                 assign_dict[seq.id] = "_".join(mnemo.generate(strength=128).split()[:6])
#             else:
#                 # calculate distances
#                 distances = []
#                 for s1, s2 in itertools.combinations(group, 2):
#                     d = pairwise_distance(
#                         s1=s1,
#                         s2=s2,
#                         vgene_field=vgene_key,
#                         jgene_field=jgene_key,
#                         cdr3_field=cdr3_key,
#                         mutations_field="mutations",
#                         shared_mutation_bonus=shared_mutation_bonus,
#                         length_penalty_multiplier=length_penalty_multiplier,
#                     )
#                     distances.append(d)
#                 # cluster
#                 linkage_matrix = fastcluster.linkage(
#                     distances,
#                     method="average",
#                     preserve_input=False,
#                 )
#                 cluster_list = fcluster(
#                     linkage_matrix,
#                     distance_cutoff,
#                     criterion="distance",
#                 )
#                 # rename clusters
#                 cluster_ids = set(cluster_list)
#                 cluster_names = {
#                     c: "_".join(mnemo.generate(strength=128).split()[:6])
#                     for c in cluster_ids
#                 }
#                 renamed_clusters = [cluster_names[c] for c in cluster_list]
#                 # assign sequences
#                 for seq, name in zip(group, renamed_clusters):
#                     assign_dict[seq.id] = name

#     # return just the assignment dict if requested
#     if return_assignment_dict:
#         return assign_dict

#     # otherwise, add the lineage name and size to the sequence DataFrame
#     lineage_size_dict = Counter(assign_dict.values())
#     lineages = [assign_dict.get(s, None) for s in df["sequence_id"]]
#     lineage_sizes = [lineage_size_dict.get(lineage, 1) for lineage in lineages]
#     df = df.with_columns(
#         pl.Series(name=lineage_field, values=lineages),
#         pl.Series(name=lineage_size_field, values=lineage_sizes),
#     )

#     # output
#     if output_path is not None:
#         if output_fmt == "airr":
#             pl.write_csv(df, output_path, separator="\t")
#         elif output_fmt == "parquet":
#             pl.write_parquet(df, output_path)
#     else:
#         return from_polars(df)


# def pairwise_distance(
#     s1: Sequence,
#     s2: Sequence,
#     shared_mutation_bonus: float = 0.65,
#     length_penalty_multiplier: Union[int, float] = 2,
#     vgene_field: str = "v_gene",
#     jgene_field: str = "j_gene",
#     cdr3_field: str = "cdr3_aa",
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

#     vgene_field : str, optional
#         Name of the field in `s1` and `s2` containing the V-gene name,
#         by default ``"v_gene"``

#     jgene_field : str, optional
#         Name of the field in `s1` and `s2` containing the J-gene name,
#         by default ``"j_gene"``

#     cdr3_field : str, optional
#         Name of the field in `s1` and `s2` containing the CDR3 sequence,
#         by default ``"cdr3_aa"``

#     mutations_field : str, optional
#         Name of the field in `s1` and `s2` containing mutation information,
#         by default ``"mutations"`` Note that this field should be a list of
#         mutation strings

#     Returns
#     -------
#     distance: float
#         The adjusted Levenshtein distance between the two sequences
#     """
#     germline_penalty = 0
#     if s1[vgene_field] != s2[vgene_field]:
#         germline_penalty += 10
#     if s1[jgene_field] != s2[jgene_field]:
#         germline_penalty += 5
#     dist = levenshtein_distance(s1[cdr3_field], s2[cdr3_field])
#     s1_len = len(s1[cdr3_field])
#     s2_len = len(s2[cdr3_field])
#     length_penalty = abs(s1_len - s2_len) * length_penalty_multiplier
#     length = min(s1_len, s2_len)
#     mutation_bonus = (
#         len(set(s1[mutations_field]) & set(s2[mutations_field])) * shared_mutation_bonus
#     )
#     score = (dist + length_penalty - mutation_bonus) / length
#     return max(score, 0.001)  # distance values can't be negative


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
