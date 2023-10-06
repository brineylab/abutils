#!/usr/bin/env python
# filename: similarity.py


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


from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd

from ..core.sequence import Sequence
from ..core.pair import Pair




def repertoire_similarity(
    repertoire1: Union[Iterable[Sequence], Iterable[Pair]],
    repertoire2: Union[Iterable[Sequence], Iterable[Pair]],
    method: str = 'morisita_horn',
    features: Union[str, Iterable[str], None] = None,
    subsample_size: Optional[int] = None,
    with_replacement: bool = False,
    n_iters=10,
    return_all_iterations: bool = False,
    pairs_only: bool = False,
    chain: Optional[str] = None,
    seed: Union[int, float] = 42,
) -> Union[float, Iterable[float]]:
    """
    
    """
     # features
    if features is None:
        features = ["v_gene", "j_gene", "cdr3_length"]
    elif isinstance(features, str):
        features = [features]
    # chain/locus
    chain_dict = {
        "heavy": ["IGH"],
        "kappa": ["IGK"],
        "lambda": ["IGL"],
        "light": ["IGK", "IGL"],
        "all": ["IGH", "IGK", "IGL"],
    }
    locus = chain_dict.get(chain.lower(), chain_dict["all"])
    # filter repertoires
    features1 = _get_features(
        repertoire=repertoire1, 
        features=features, 
        locus=locus, 
        pairs_only=pairs_only
        )
    features2 = _get_features(
        repertoire=repertoire2, 
        features=features, 
        locus=locus, 
        pairs_only=pairs_only
        )
    # calculate similarity
    np.random.seed(seed)
    similarities = []
    similarity_funcs = {
        'morisita_horn': morisita_horn,
        'kullback-leibler': kullback_leibler,
        'jensen-shannon': jensen_shannon,
        'jaccard': jaccard_similarity,
        'bray-curtis': bray_curtis,
        'renkonen': renkonen,
        'cosine': cosine_similarity,
    }
    similarity_func = similarity_funcs.get(method.lower(), "morisita_horn")
    for i in range(n_iters):
        # subsample
        if subsample_size is None:
            subsample_size = min(len(features1), len(features2))
        subsample1 = np.random.choice(features1, size=subsample_size, replace=with_replacement)
        subsample2 = np.random.choice(features2, size=subsample_size, replace=with_replacement)
        # count features
        counts1 = pd.Series(subsample1).value_counts()
        counts2 = pd.Series(subsample2).value_counts()
        df = pd.DataFrame({'sample1': counts1, 'sample2': counts2}).fillna(0)
        # similarity
        similarity = similarity_func(df['sample1'], df['sample2'])
        similarities.append(similarity)
    if return_all_iterations:
        return similarities
    else:
        return np.mean(similarities)





def _get_features(
        repertoire: Union[Iterable[Sequence], Iterable[Pair]],
        features: Iterable[str],
        locus: Iterable[str],
        pairs_only: bool = False,
) -> Iterable[str]:
    """
    
    """
    features = []
    # process Pair objects
    if all([isinstance(r, Pair) for r in repertoire]):
        if pairs_only:
            repertoire = [r for r in repertoire if r.is_paired]
        for r in repertoire:
            feat_list = []
            if "IGH" in locus and r.heavy is not None:
                for f in features:
                    feat = r.heavy["f"]
                    feat_list.append(feat if feat is not None else "NA")
            if r.light is not None and r.light["locus"] in locus:
                for f in features:
                    feat = r.light["f"]
                    feat_list.append(feat if feat is not None else "NA")
            features.append("__".join(feat_list))
    else:
        for r in repertoire:
            feat_list = []
            if r["locus"] in locus:
                for f in features:
                    feat = r[f]
                    feat_list.append(feat if feat is not None else "NA")
            features.append("__".join(feat_list))
    return features
            



#-----------------------------------
#       SIMILARITY METHODS
#-----------------------------------

def morisita_horn(sample1, sample2):
    """
    
    """
    X = sum(sample1)
    Y = sum(sample2)
    XY = X * Y
    sumXiYi = 0
    sumXiSq = 0
    sumYiSq = 0
    for x, y in zip(sample1, sample2):
        sumXiYi += x * y
        sumXiSq += x * x
        sumYiSq += y * y
    numerator = 2 * sumXiYi
    denominator = (float(sumXiSq) / (X * X) + float(sumYiSq) / (Y * Y)) * XY
    return numerator / denominator





