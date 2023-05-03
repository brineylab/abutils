from collections import Counter
import itertools
import sys
from typing import Iterable, Union

import pandas as pd
import numpy as np

from Levenshtein import distance

import fastcluster as fc
from scipy.cluster.hierarchy import fcluster

from mnemonic import Mnemonic

from ..core.sequence import Sequence
from ..tools import cluster
from ..utils.utilities import nested_dict_lookup


def clonify(
    sequences: Iterable[Sequence],
    distance_cutoff: float = 0.32,
    shared_mutation_bonus: float = 0.65,
    length_penalty_multiplier: Union[int, float] = 2,
    group_by_v: bool = True,
    group_by_j: bool = True,
    preclustering: bool = False,
    preclustering_threshold: float = 0.65,
    preclustering_field: str = "cdr3_nt",
    lineage_field: str = "lineage",
    lineage_size_field: str = "lineage_size",
    annotation_format: str = "airr",
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
    sequences : Iterable[Sequence]
        Iterable of ``abutils.core.Sequence`` objects containing annotated sequence data.  

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
        ``group_by_j == True``, sequences will be grouped by both V-gene and J-gene use.
        
    group_by_j : bool, default=True
        If ``True``, sequences are grouped by J-gene use prior to lineage assignment. This 
        option is additive with ``group_by_v``. For example, if ``group_by_v == True`` and
        ``group_by_j == True``, sequences will be grouped by both V-gene and J-gene use.

    preclustering : bool, default=False  
        If ``True``, V/J groups are pre-clustered on the `preclustering_field` sequence, 
        which can potentially speed up lineage assignment and reduce memory usage. 
        If ``False``, each V/J group is processed in its entirety without pre-clustering. 

    preclustering_threshold : float, default=0.65  
        Identity threshold for pre-clustering the V/J groups prior to lineage assignment. 

    preclustering_field : str, default='cdr3_nt'  
        Annotation field on which to pre-cluster sequences.  

    lineage_field : str, default='lineage'  
        Name of the lineage assignment field.  

    lineage_size_field : str, default='lineage_size'  
        Name of the lineage size field.  

    annotation_format : str, default='airr'  
        Format of the input sequence annotations. Choices are ``'airr'`` or ``'json'``.  
        
    return_assignment_dict : bool, default=False  
        If ``True``, a dictionary linking sequence IDs to lineage names will be returned. 
        If ``False``, the input ``anndata.AnnData`` object will be returned, with lineage 
        annotations included.  

    Returns
    -------
    output : ``Iterable[Sequence]`` or ``dict``
        By default (``return_assignment_dict == False``), a list of input sequences are \
        returned with two additional annotations - the lineage name and lineage size. \
        If ``return_assignment_dict == True``, \
        a ``dict`` mapping sequence IDs to lineage names is returned. 


    .. _clonify: https://github.com/briney/clonify    
    """
    # select the appropriate data fields
    if annotation_format.lower() == "airr":
        vgene_key = "v_gene"
        jgene_key = "j_gene"
        cdr3_key = "cdr3_aa"
        muts_key = "v_mutations"
    elif annotation_format.lower() == "json":
        vgene_key = "v_gene.gene"
        jgene_key = "j_gene.gene"
        cdr3_key = "cdr3_aa"
        muts_key = "var_muts_nt.muts"
    else:
        error = "ERROR: "
        error += f'annotation_format must be either "airr" or "json", but you provided {annotation_format}'
        print("\n")
        print(error)
        print("\n")
        sys.exit()
    # group sequences by V/J genes
    vj_group_dict = {}
    heavies = [s for s in sequences if s["locus"] == "IGH"]
    for h in heavies:
        # build new Sequence objects using just the data we need
        s = Sequence(h.sequence, id=h.id)
        s["v_gene"] = nested_dict_lookup(h, vgene_key.split("."))
        s["j_gene"] = nested_dict_lookup(h, jgene_key.split("."))
        s["cdr3"] = nested_dict_lookup(h, cdr3_key.split("."))
        if annotation_format.lower() == "json":
            muts = nested_dict_lookup(h, muts_key.split("."), [])
            s["mutations"] = [f"{m['position']}:{m['was']}>{m['is']}" for m in muts]
        else:
            muts = h[muts_key]
            if any([pd.isnull(muts), muts is None]):
                s["mutations"] = []
            else:
                muts = muts.split("|")
                s["mutations"] = [m for m in muts if m.strip()]
        required_fields = ["v_gene", "j_gene", "cdr3", "mutations"]
        if preclustering:
            s["preclustering"] = nested_dict_lookup(h, preclustering_field.split("."))
            required_fields.append("preclustering")
        if any([s[v] is None for v in required_fields]):
            continue
        # group sequences by VJ gene use
        # vj = f"{s['v_gene']}__{s['j_gene']}"
        vj = ""
        if group_by_v:
            vj += s["v_gene"]
        if group_by_j:
            vj += s["j_gene"]
        if vj not in vj_group_dict:
            vj_group_dict[vj] = []
        vj_group_dict[vj].append(s)
    # assign lineages
    assignment_dict = {}
    mnemo = Mnemonic("english")
    for vj_group in vj_group_dict.values():
        # preclustering
        if preclustering:
            seq_dict = {s.id: s for s in vj_group}
            cluster_seqs = [Sequence(s[preclustering_field], id=s.id) for s in vj_group]
            clusters = cluster(cluster_seqs, threshold=preclustering_threshold)
            groups = [[seq_dict[i] for i in c.seq_ids] for c in clusters]
        else:
            groups = [
                vj_group,
            ]
        for group in groups:
            if len(group) == 1:
                seq = group[0]
                assignment_dict[seq.id] = "_".join(
                    mnemo.generate(strength=128).split()[:6]
                )
                continue
            # build a distance matrix
            dist_matrix = []
            for s1, s2 in itertools.combinations(group, 2):
                d = _clonify_distance(
                    s1, s2, shared_mutation_bonus, length_penalty_multiplier
                )
                dist_matrix.append(d)
            # cluster
            linkage_matrix = fc.linkage(
                dist_matrix, method="average", preserve_input=False
            )
            cluster_list = fcluster(
                linkage_matrix, distance_cutoff, criterion="distance"
            )
            # rename clusters
            cluster_ids = list(set(cluster_list))
            cluster_names = {
                c: "_".join(mnemo.generate(strength=128).split()[:6])
                for c in cluster_ids
            }
            renamed_clusters = [cluster_names[c] for c in cluster_list]
            # assign sequences
            for seq, name in zip(vj_group, renamed_clusters):
                assignment_dict[seq.id] = name
    lineage_size_dict = Counter(assignment_dict.values())
    # return assignments
    if return_assignment_dict:
        return assignment_dict
    for s in input_seqs:
        if l := assignment_dict.get(s.id, None) is not None:
            s[lineage_field] = l
            s[lineage_size_field] = lineage_size_dict.get(l, 1)
    return input_seqs
    # lineage_assignments = [assignment_dict.get(n, np.nan) for n in adata.obs_names]
    # lineage_sizes = [lineage_size_dict.get(l, np.nan) for l in lineage_assignments]
    # adata.obs[lineage_field] = lineage_assignments
    # adata.obs[lineage_size_field] = lineage_sizes
    # return adata


def pairwise_distance(
    s1: Sequence,
    s2: Sequence,
    shared_mutation_bonus: float = 0.65,
    length_penalty_multiplier: Union[int, float] = 2,
    cdr3_field: str = "cdr3",
    mutations_field: str = "mutations",
) -> float:
    """Computes length and mutation adjusted Levenshtein distance for a pair of sequences.

    Parameters
    ----------
    s1 : abutils.Sequence
        input sequence
    s2 : abutils.Sequence
        input sequence
    shared_mutation_bonus : float, optional
        The bonus for each shared mutation, by default 0.65
    length_penalty_multiplier : Union[int, float], optional
        Used to compute the penalty for differences in CDR3 length. The length
        difference is multiplied by `length_penalty_multiplier`, by default 2
    cdr3_field : str, optional
        Name of the field in `s1` and `s2` containing the CDR3 sequence,
        by default "cdr3"
    mutations_field : str, optional
        Name of the field in `s1` and `s2` containing mutation information,
        by default "mutations"

    Returns
    -------
    float
        distance
    """
    if len(s1[cdr3_field]) == len(s2[cdr3_field]):
        dist = sum([i != j for i, j in zip(s1[cdr3_field], s2[cdr3_field])])
    else:
        dist = distance(s1[cdr3_field], s2[cdr3_field])
    length_penalty = (
        abs(len(s1[cdr3_field]) - len(s2[cdr3_field])) * length_penalty_multiplier
    )
    length = min(len(s1[cdr3_field]), len(s2[cdr3_field]))
    shared_mutations = list(set(s1[mutations_field]) & set(s2[mutations_field]))
    mutation_bonus = len(shared_mutations) * shared_mutation_bonus
    score = (dist + length_penalty - mutation_bonus) / length
    return max(score, 0.001)  # distance values can't be negative
