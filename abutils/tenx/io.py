#!/usr/bin/env python
# filename: io.py


#
# Copyright (c) 2020 Bryan Briney
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


import re

import numpy as np

import scanpy as sc


def read_10x_mtx(mtx_file, gex_only=False, hash_regex='cell ?hash', ignore_hash_regex_case=True,
                 transform_hashes=True, ignore_zero_median_hashes=True,
                 transform_features=True, ignore_zero_median_features=True):
    '''
    Reads a 10x Genomics matrix file and outputs gex, cell hash, and feature barcode data.

    Args:

        mtx_file (str): path to the 10x Genomics matrix folder (as accepted by ``scanpy.read_10x_mtx()``)

        gex_only (bool): If ``True``, return only gene expression data and ignore features and hashes.
            Default is ``False``.

        hash_regex (str): A regular expression (regex) string used to identify cell hashes. The regex 
            must be found in all hash names. The default is ``'cell ?hash'``, which combined with the
            default setting for ``ignore_hash_regex_case``, will match ``'cellhash'`` or ``'cell hash'``
            in any combination of upper and lower case letters.

        ignore_hash_regex_case (bool): If ``True``, searching for ``hash_regex`` will ignore case.
            Default is ``True``.

        transform_hashes (bool): If ``True``, cell hash counts will be log2 transformed 
            (after adding 1 to the raw count). Default is ``True``.
        
        transform_features (bool): If ``True``, feature counts will be log2 transformed 
            (after adding 1 to the raw count). Default is ``True``.

        ignore_zero_median_hashes (bool): If ``True``, any hashes containing a meadian
            count of ``0`` will be ignored and not returned in the hash dataframe. Default
            is ``True``.
        
        ignore_zero_median_features (bool): If ``True``, any features containing a meadian
            count of ``0`` will be ignored and not returned in the feature dataframe. Default
            is ``True``.

    Returns:

        gex (anndata.AnnData): An ``AnnData`` object containing gene expression data. If ``gex_only``
            is ``True``, only ``gex`` is returned.

        hash_df (pandas.DataFrame): A ``DataFrame`` containing cell hash data. Not returned if 
            ``gex_only`` is ``True``.

        feature_df (pandas.DataFrame): A ``DataFrame`` containing feature data. Not returned if 
            ``gex_only`` is ``True``.
    '''
    # read input file
    adata = sc.read_10x_mtx(mtx_file, gex_only=False)
    # split input files
    gex = adata[:,adata.var.feature_types == 'Gene Expression']
    if gex_only:
        return gex
    non_gex = adata[:,adata.var.feature_types != 'Gene Expression']
    # parse out features and hashes
    if ignore_hash_regex_case:
        pattern = re.compile(hash_regex, flags=re.IGNORECASE)
    else:
        pattern = re.compile(hash_regex)
    hashes = non_gex[:, [re.search(pattern, i) is not None for i in non_gex.var.gene_ids]]
    features = non_gex[:, [re.search(pattern, i) is None for i in non_gex.var.gene_ids]]
    # make hash DataFrame
    hash_df = sc.get.obs_df(hashes, hashes.var_names)
    if ignore_zero_median_hashes:
        hash_df = hash_df[[h for h in hash_df.columns.values if hash_df[h].median() > 0]]
    if transform_hashes:
        hash_df += 1
        hash_df = hash_df.apply(np.log2)
    # make feature DataFrame
    feature_df = sc.get.obs_df(features, features.var_names)
    if ignore_zero_median_features:
        feature_df = feature_df[[h for h in feature_df.columns.values if feature_df[h].median() > 0]]
    if transform_features:
        feature_df += 1
        feature_df = feature_df.apply(np.log2)
    return gex, hash_df, feature_df