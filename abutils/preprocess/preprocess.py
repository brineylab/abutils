#!/usr/bin/env python
# filename: preprocess.py


#
# Copyright (c) 2025 Bryan Briney & Benjamin Nemoz
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
import time
import datetime
from tqdm.auto import tqdm

import polars as pl
from ..io import list_files, make_dir, to_fasta, from_polars
from ..tools import cluster, alignment


__all__ = ['deduplicate', 'reduction']


def deduplicate(project_folder, 
                  output=None, 
                  pool=True, 
                  umi=False, 
                  keep_read_numbers=False,
                  read_number_separator='|',
                  large_files=False, 
                  debug=False):
    '''
    A polyvalent tool for deduplication of assigned reads. This function takes as input the AIRR-compliant 
    tables and is specifically designed to handle extremely large files with a minimal footprint.

    Parameters:
    project_folder (str): Path to the project folder containing AIRR-compliant tables.
    output (str, optional): Subdirectory for output files. Created if non-existent. Defaults to None.
    pool (bool, optional): If True, pool all samples together. Defaults to True.
    umi (bool, optional): If True, use UMI for deduplication. Defaults to False.
    keep_read_numbers (bool, optional): If True, read numbers will be added to sequence names. Defaults to False.
    large_files (bool, optional): If True, optimize for large files (>100Go). Defaults to False.
    debug (bool, optional): If True, print debug information. Defaults to False.

    Returns: None. A fasta file is written to disk.
    '''

    start = time.time()

    if os.path.isfile(project_folder):
        project_folder = os.path.dirname(project_folder)
        if debug:
            print(f"Project folder is set to: {project_folder}")
         
    files = sorted([f for f in list_files(project_folder, extension='tsv', recursive=True, ) if not 'tmp' in f])

    if len(files) == 0:
        print('No files found. Exiting.')
    elif len(files) == 1:
        pool = False

    if output:
        project_folder = os.path.join(project_folder, output)
        make_dir(project_folder)
    

    total_sequences = 0
    pooled = []
    
    for file in files:
        sample = os.path.basename(file).split('.')[0]

        print(f"Processing {sample}")
        print("-"*(len(sample)+11)+'\n')
        
        if umi:
            df = pl.read_csv(file, columns=['sequence_id', 'sequence', 'umi'], separator='\t', null_values="None", low_memory=True if large_files else False, )
        else:
            df = pl.read_csv(file, columns=['sequence_id', 'sequence'], separator='\t', null_values="None", low_memory=True if large_files else False, )

        print(f"Loaded {df.height:,} annotated sequences")
        total_sequences += df.shape[0]

        if umi:
            if debug:
                print(f"Deduplicating with UMI ({df.unique(subset=['umi']).height} unique UMIs)")
            
            if not keep_read_numbers:
                df_unique = df.unique(subset=["sequence", "umi"]).sort("sequence")
            else:
                df_unique = (
                                df.group_by(["sequence", "umi"]).agg(
                                    pl.count().alias("count"), 
                                    pl.col("sequence_id").first(), 
                                )
                                .sort("sequence")
                            )
        else:
            if not keep_read_numbers:
                df_unique = df.unique(subset=["sequence"]).sort("sequence")
            else:
                df_unique = (
                                df.group_by("sequence").agg(
                                    pl.count().alias("count"), 
                                    pl.col("sequence_id").first(), 
                                )
                                .sort("sequence")
                            )

        msg = f"Found {df_unique.height:,} unique sequences"
        if pool:
            msg += " added to the pool\n"
        print(msg)

        if pool:
            pooled.append(df_unique)
        else:
            fasta_file = os.path.join(project_folder, sample)+'.fasta'
            if keep_read_numbers:
                df_unique = df_unique.with_columns(
                                                    (pl.col("sequence_id") + read_number_separator + pl.col("count").cast(pl.Utf8)).alias("sequence_id")
                                                )

            tuples_list = list(zip(df_unique["sequence_id"].to_list(), df_unique["sequence"].to_list()))
            to_fasta(tuples_list, fasta_file, )
            print(f"Output written to fasta file: {fasta_file}\n")

    if pool:
        pool_df = pl.concat(pooled)

        if not keep_read_numbers:
            pool_unique = pool_df.unique(subset=["sequence"]).sort("sequence")
        else:
            pool_unique = (
                                pool_df.group_by("sequence").agg(
                                    pl.count().alias("count"), 
                                    pl.col("sequence_id").first(), 
                                )
                                .sort("sequence")
                            )
                            
        print(f"\nFound {pool_unique.height:,} unique sequences in pooled data")
    
        fasta_file = os.path.join(project_folder, "deduplicated_pool.fasta")
        
        tuples_list = list(zip(pool_unique["sequence_id"].to_list(), pool_unique["sequence"].to_list()))

        to_fasta(tuples_list, fasta_file, )
        print(f"Output written to fasta file: {fasta_file}\n")

    end = time.time()
    duration = end - start
    print(f"Deduplication complete. Total elapsed time: {datetime.timedelta(seconds=int(duration))}.")

    return



def reduction(project_folder, 
              output=None, 
              pool=True, 
              umi=False, 
              keep_cluster_sizes=False, 
              cluster_sizes_separator='|',
              min_cluster_size=3, 
              clustering_threshold=0.975, 
              consentroid='centroid', 
              large_files=False, 
              debug=False):
    """
    This function takes as an input AIRR-compliant tables (tsv) and proceeds to
    data reduction by clustering sequences to a high identity threshold.

    This is specifically designed to handle large files with minimal footprint.
    Preclustering can be applied to increase performance over large datasets.

    Parameters:
    project_folder (str): Path to the project folder containing AIRR-compliant tables.
    output (str, optional): Subdirectory for output files. Created if non-existent. Defaults to None.
    pool (bool, optional): If True, pool all samples together. Defaults to True.
    umi (bool, optional): If True, use UMI for clustering. Defaults to False.
    keep_cluster_sizes (bool, optional): If True, cluster sizes will be added to sequence names. Defaults to False.
    cluster_sizes_separator (str, optional): Separator for cluster sizes in sequence names. Defaults to '|'.
    min_cluster_size (int, optional): Minimum cluster size to consider. Defaults to 3.
    clustering_threshold (float, optional): Identity threshold for clustering. Defaults to 0.975.
    consentroid (str, optional): Method to determine cluster representative ('centroid' or 'consensus'). Defaults to 'centroid'.
    large_files (bool, optional): If True, optimize for large files (>100Go). Defaults to False.
    debug (bool, optional): If True, print debug information. Defaults to False.

    Returns: None. A fasta file is written to disk.
    """

    start = time.time()

    if output:
        project_folder = os.path.join(project_folder, output)
        make_dir(project_folder)

    files = sorted([f for f in list_files(project_folder, extension='tsv', recursive=True, ) if 'airr' in f])

    total_sequences = 0

    if pool:
        pooled_heavies = []
        pooled_lights = []

    for file in files:
        sample = os.path.basename(file).split('.')[0]

        print(f"Processing {sample}")
        print("-"*(len(sample)+11)+'\n')
        
        keys = ['sequence_id', 'v_gene', 'j_gene', 'locus', 'sequence']
        if umi:
            keys.append('umi')

        df = pl.read_csv(file, columns=keys, separator='\t', null_values="None", low_memory=True if large_files else False, )

        print(f"Loaded {df.height:,} annotated sequences")
        total_sequences += df.height

        heavies = df.filter(pl.col('locus') == 'IGH')
        lights = df.filter(pl.col('locus') != 'IGH')

        if pool:
            pooled_heavies.append(heavies)
            pooled_lights.append(lights)
            
        else:
            sample_consentroids = []

            # Process heavy chains
            heavies = heavies.with_columns((pl.col("v_gene") + "_" + pl.col("j_gene")).alias("vj_bin"))
            heavy_consentroids = []

            for vj in tqdm(heavies['vj_bin'].unique().to_list(), desc='Heavy chains'):
                group = heavies.filter(pl.col('vj_bin') == vj)
                
                if group.height < min_cluster_size:
                    continue

                group_sequences = from_polars(group)

                clusters = cluster.cluster(group_sequences, threshold=clustering_threshold)
                for c in clusters:
                    if c.size < min_cluster_size:
                        continue
                    else:
                        if consentroid == 'centroid':
                            heavy_consentroids.append(c.centroid)
                        elif consentroid == 'consensus':
                            _consensus = alignment.make_consensus([s for s in c.sequences], algo='mafft', alignment_kwargs={'mafft_bin':'mafft','threads':1}, debug=debug)
                            heavy_consentroids.append(_consensus)

            sample_consentroids.extend(heavy_consentroids)

            # Process light chains
            lights = lights.with_columns((pl.col("v_gene") + "_" + pl.col("j_gene")).alias("vj_bin"))
            light_consentroids = []

            for vj in tqdm(lights['vj_bin'].unique().to_list(), desc='Light chains'):
                group = lights.filter(pl.col('vj_bin') == vj)
                
                if group.height < min_cluster_size:
                    continue

                group_sequences = from_polars(group)

                clusters = cluster.cluster(group_sequences, threshold=clustering_threshold)
                for c in clusters:
                    if c.size < min_cluster_size:
                        continue
                    else:
                        if consentroid == 'centroid':
                            light_consentroids.append(c.centroid)
                        elif consentroid == 'consensus':
                            _consensus = alignment.make_consensus([s for s in c.sequences], algo='mafft', alignment_kwargs={'mafft_bin':'mafft','threads':1}, debug=debug)
                            light_consentroids.append(_consensus)

            sample_consentroids.extend(light_consentroids)

            fasta_file = os.path.join(project_folder, sample)+'_reduced.fasta'
            to_fasta(sample_consentroids, fasta_file)

    if pool:
        heavies = pl.concat(pooled_heavies)
        lights = pl.concat(pooled_lights)

        all_consentroids = []

        # Process heavy chains
        heavies = heavies.with_columns((pl.col("v_gene") + "_" + pl.col("j_gene")).alias("vj_bin"))
        heavy_consentroids = []

        for vj in tqdm(heavies['vj_bin'].unique().to_list(), desc='Heavy chains'):
            group = heavies.filter(pl.col('vj_bin') == vj)
            
            if group.height < min_cluster_size:
                continue

            group_sequences = from_polars(group)

            clusters = cluster.cluster(group_sequences, threshold=clustering_threshold)
            for c in clusters:
                if c.size < min_cluster_size:
                    continue
                else:
                    if consentroid == 'centroid':
                        heavy_consentroids.append(c.centroid)
                    elif consentroid == 'consensus':
                        _consensus = alignment.make_consensus([s for s in c.sequences], algo='mafft', alignment_kwargs={'mafft_bin':'mafft','threads':1}, debug=debug)
                        heavy_consentroids.append(_consensus)

        all_consentroids.extend(heavy_consentroids)

        # Process light chains
        lights = lights.with_columns((pl.col("v_gene") + "_" + pl.col("j_gene")).alias("vj_bin"))
        light_consentroids = []

        for vj in tqdm(lights['vj_bin'].unique().to_list(), desc='Light chains'):
            group = lights.filter(pl.col('vj_bin') == vj)
            
            if group.height < min_cluster_size:
                continue

            group_sequences = from_polars(group)

            clusters = cluster.cluster(group_sequences, threshold=clustering_threshold)
            for c in clusters:
                if c.size < min_cluster_size:
                    continue
                else:
                    if consentroid == 'centroid':
                        light_consentroids.append(c.centroid)
                    elif consentroid == 'consensus':
                        _consensus = alignment.make_consensus([s for s in c.sequences], algo='mafft', alignment_kwargs={'mafft_bin':'mafft','threads':1}, debug=debug)
                        light_consentroids.append(_consensus)

        all_consentroids.extend(light_consentroids)

        fasta_file = os.path.join(project_folder, "reduced_pool.fasta")
        to_fasta(all_consentroids, fasta_file)

    end = time.time()
    duration = end - start
    print(f"Data reduction complete. Total elapsed time: {datetime.timedelta(seconds=int(duration))}.")

    return 