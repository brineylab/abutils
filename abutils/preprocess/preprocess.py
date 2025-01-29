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

import polars as pl
from .io import list_files, make_dir, to_fasta


__all__ = ['deduplication', 'reduction']


def deduplication(project_folder, 
                  output=None, 
                  pool=True, 
                  umi=False, 
                  keep_read_numbers=False, 
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

    Returns:
    None
    '''

    start = time.time()
    files = sorted([f for f in list_files(project_folder, extension='tsv', recursive=True, ) if 'airr' in f])

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

        print(f"Loaded {df.shape[0]:,} annotated sequences")
        total_sequences += df.shape[0]

        df_unique = df.unique(subset=["sequence"]).sort("sequence")
        # TO-DO: deduplication with UMI
        # TO-DO: count-persistent deduplication

        msg = f"Found {df_unique.shape[0]:,} unique sequences"
        if pool:
            msg += " added to the pool\n"
        print(msg)

        if pool:
            pooled.append(df_unique)
        else:
            fasta_file = os.path.join(project_folder, sample)+'.fasta'
            tuples_list = list(zip(df["sequence_id"].to_list(), df["sequence"].to_list()))
            to_fasta(tuples_list, fasta_file, )
            print(f"Output written to fasta file: {fasta_file}\n")

    if pool:
        pool_df = pl.concat(pooled)
        pool_unique = pool_df.unique(subset=["sequence"]).sort("sequence")
        
        print(f"\nFound {pool_unique.shape[0]:,} unique sequences in pooled data")
    
        fasta_file = os.path.join(project_folder, "deduplicated_pool.fasta")
        tuples_list = list(zip(pool_unique["sequence_id"].to_list(), pool_unique["sequence"].to_list()))
        to_fasta(tuples_list, fasta_file, )
        print(f"Output written to fasta file: {fasta_file}\n")

    end = time.time()
    duration = end - start
    print(f"Deduplication complete. Total elapsed time: {datetime.timedelta(seconds=duration)}.")

    return



def reduction(project_folder, 
              output=None, 
              pool=True, 
              umi=False, 
              preclustering=False,
              keep_cluster_sizes=False, 
              min_cluser_size=3, 
              clustering_threshold=0.975, 
              consentroid='centroid', 
              large_files=False, 
              debug=False):
    """
    This function takes as an input AIRR-compliant tables (tsv) and proceed to
    data reduction by clustering sequences to a high identity threshold

    This is specifically designed to handle large files with minimal footprint
    Preclustering can be applied to increase performance over large datasets
    
    Returns a fasta file of consensus or centroid sequences
    """

    start = time.time()
    files = sorted([f for f in list_files(project_folder, extension='tsv', recursive=True, ) if 'airr' in f])

    # TO-DO: everything!

    end = time.time()
    duration = end - start
    print(f"Data reduction complete. Total elapsed time: {datetime.timedelta(seconds=duration)}.")


    return 