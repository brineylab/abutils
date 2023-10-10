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


import itertools
import math
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd

from ..core.sequence import Sequence
from ..core.pair import Pair


__all__ = [
    "RepertoireSimilarity",
    "RepertoireSimilarities",
    "repertoire_similarity",
    # "normalize",
    # "make_continuous",
    "morisita_horn",
    "kullback_leibler",
    "jensen_shannon",
    "jaccard_similarity",
    "renkonen",
    "bray_curtis",
    "cosine_similarity",
]


class RepertoireSimilarity:
    """
    Object containing repertoire similarity data for a single pairwise comparison.
    """

    def __init__(
        self,
        similarities: Iterable[float],
        method: str,
        sample1_name: Optional[str] = None,
        sample2_name: Optional[str] = None,
    ):
        self.similarities = similarities
        self.method = method
        self.sample1_name = sample1_name if sample1_name is not None else "sample1"
        self.sample2_name = sample2_name if sample2_name is not None else "sample2"
        self._df = None

    def __repr__(self):
        return "<RepertoireSimilarity: {}>".format(self.method)

    def __str__(self):
        return "<RepertoireSimilarity: {}>".format(self.method)

    def __len__(self):
        return len(self.similarities)

    @property
    def df(self):
        """
        Returns a pandas ``DataFrame`` of similarity data.
        """
        if self._df is None:
            data = []
            for i, s in enumerate(self.similarities, 1):
                d = {
                    "sample1": self.sample1_name,
                    "sample2": self.sample2_name,
                    "method": self.method,
                    "iteration": i,
                    "similarity": s,
                }
                data.append(d)
            self._df = pd.DataFrame(data)
        return self._df

    @property
    def mean(self):
        """
        Returns the mean similarity value (if multiple iterations were performed)
        """
        return np.mean(self.similarities)

    @property
    def median(self):
        """
        Returns the median similarity value (if multiple iterations were performed)
        """
        return np.median(self.similarities)

    @property
    def min(self):
        """
        Returns the minimum similarity value (if multiple iterations were performed)
        """
        return np.min(self.similarities)

    @property
    def max(self):
        """
        Returns the maximum similarity value (if multiple iterations were performed)
        """
        return np.max(self.similarities)

    @property
    def std(self):
        """
        Returns the similarity standard deviation (if multiple iterations were performed)
        """
        return np.std(self.similarities)

    @property
    def sem(self):
        """
        Returns the similarity standard error of means (if multiple iterations were performed)
        """
        return np.std(self.similarities) / np.sqrt(len(self.similarities))


class RepertoireSimilarities:
    """
    Object containing repertoire similarity data for a multiple pairwise comparisons.
    """

    def __init__(self, similarities: Optional[Iterable[RepertoireSimilarity]] = None):
        self.similarities = similarities if similarities is not None else []
        self.df_needs_update = True
        self._df = None

    @property
    def df(self):
        """
        Returns a pandas ``DataFrame`` of similarity data.
        """
        if self._df is None or self.df_needs_update:
            if self.similarities:
                self._df = pd.concat([s.df for s in self.similarities])
                self.df_needs_update = False
        return self._df

    @property
    def means(self):
        """
        Returns an array of mean similarity values, one for each pairwise comparison.
        """
        return np.array([s.mean for s in self.similarities])

    @property
    def medians(self):
        """
        Returns an array of median similarity values, one for each pairwise comparison.
        """
        return np.array([s.median for s in self.similarities])

    @property
    def mins(self):
        """
        Returns an array of minimum similarity values, one for each pairwise comparison.
        """
        return np.array([s.min for s in self.similarities])

    @property
    def maxs(self):
        """
        Returns an array of maximum similarity values, one for each pairwise comparison.
        """
        return np.array([s.max for s in self.similarities])

    @property
    def stds(self):
        """
        Returns an array of standard deviations, one for each pairwise comparison.
        """
        return np.array([s.std for s in self.similarities])

    @property
    def sems(self):
        """
        Returns an array of standard errors of the mean, one for each pairwise comparison.
        """
        return np.array([s.sem for s in self.similarities])

    def squareform(self, method=None, agg=np.mean):
        """
        Returns a squareform representation of the similarity data.

        Parameters
        ----------
        method : str, optional
            The similarity method to use, if multiple similarity methods are present in the
            ``RepertoireSimilarities`` object. If ``None``, and multiple methods exist, all methods
            will be pooled (not recommended). Default is ``None``, which is the preferred
            argument when the data contains only a single similarity method.

        agg : function, optional
            The function used to aggregate similarity values. Default is ``np.mean``.
        """
        sq_data = {}
        for s1 in self.df["sample1"].unique():
            for s2 in self.df["sample2"].unique():
                _df = self.df[(self.df["sample1"] == s1) & (self.df["sample2"] == s2)]
                if method is not None:
                    _df = _df[_df["method"] == method]
                if _df.shape[0] == 0:
                    continue
                similarity = agg(_df["similarity"])
                if s1 not in sq_data:
                    sq_data[s1] = {}
                sq_data[s1][s2] = similarity
        return pd.DataFrame(sq_data)

    def add(self, similarity):
        """
        Adds a RepertoireSimilarity object.

        Parameters
        ----------
        similarity : RepertoireSimilarity
            A RepertoireSimilarity object.

        """
        self.similarities.append(similarity)
        self.df_needs_update = True


def repertoire_similarity(
    repertoires: Iterable[Union[Iterable[Sequence], Iterable[Pair]]],
    names: Optional[Iterable[str]] = None,
    method: str = "morisita-horn",
    features: Union[str, Iterable[str], None] = None,
    n_iters: int = 1,
    subsample_size: Optional[int] = None,
    sample_with_replacement: bool = False,
    pairs_only: bool = False,
    chain: Optional[str] = None,
    force_self_comparisons: bool = False,
    seed: Union[int, float, str, None] = None,
) -> Union[float, RepertoireSimilarity, RepertoireSimilarities]:
    """
    Compute the pairwise similarity between two or more repertoires.

    Parameters
    ----------
    repertoires : list
        A list of repertoires. Each repertoire can be a list of ``abutils.Sequence``
        objects or ``abutils.Pair`` objects. Individual repertoires cannot mix ``Sequence``
        and ``Pair`` objects, but the list of repertoires can contain both.

    names : list, optional
        A list of names for the repertoires. If ``None``, names will be generated as
        ``"repertoire1"``, ``"repertoire2"``, etc.

    method : str, optional
        The similarity method to use. Default is ``"morisita-horn"``.

    features : str, list, optional
        The features to use for similarity calculation. Default is ``["v_gene", "j_gene", "cdr3_length"]``.

    n_iters : int, optional
        The number of iterations to perform. Default is ``1``.

    subsample_size : int, optional
        The number of sequences to subsample from each repertoire. If ``None``, the smallest repertoire
        will be used as the subsample size. Default is ``None``.

    sample_with_replacement : bool, optional
        Whether to subsample with replacement. Default is ``False``.

    pairs_only : bool, optional
        Whether to use only paired sequences. Default is ``False``.

    chain : str, optional
        The chain to use for similarity calculation. Default is ``None``, which will use all chains.

    force_self_comparisons : bool, optional
        Whether to force self-comparisons. Only used if exactly two repertoires are provided.
        Default is ``False``, which performs only a single pairwise comparison when exactly
        two repertoires are provided. If more than two repertoires are provided, all pairwise
        comparisons will be performed, including self-comparisons.

    seed : int, float, str, optional
        The seed to use for random number generation. Default is ``None``.

    Returns
    -------
    float or RepertoireSimilarity or RepertoireSimilarities
        If only two repertoires are provided and `n_iters` is ``1``, the similarity value will be returned as a float.
        If only two repertoires are provided and `n_iters` is greater than ``1``, a ``RepertoireSimilarity`` object will be
        returned. If more than two repertoires are provided, a ``RepertoireSimilarities`` object will be
        returned.

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
    if chain is None:
        locus = chain_dict["all"]
    elif isinstance(chain, str):
        locus = chain_dict.get(chain.lower(), [chain])
    else:
        locus = chain
    # subsample size
    if subsample_size is None:
        subsample_size = min([len(r) for r in repertoires])
    # names
    if names is None:
        names = [f"repertoire{i}" for i, _ in enumerate(repertoires, 1)]
    elif len(names) != len(repertoires):
        err = f"\nThe number of names must match the number of repertoires.\n"
        err += f"You provided {len(names)} names and {len(repertoires)} repertoires.\n"
        raise RuntimeError(err)
    # process repertoires
    if (n := len(repertoires)) < 2:
        err = f"\nYou must provide at least two repertoires for comparison. You provided {n}\n"
        raise RuntimeError(err)
    rep_dict = {n: r for n, r in zip(names, repertoires)}
    # similarities
    similarities = []
    if len(names) == 2 and not force_self_comparisons:
        comps = [names]
    else:
        comps = itertools.product(names, repeat=2)
    for name1, name2 in comps:
        repertoire1 = rep_dict[name1]
        repertoire2 = rep_dict[name2]
        sim = _calculate_similarity(
            repertoire1=repertoire1,
            repertoire2=repertoire2,
            name1=name1,
            name2=name2,
            method=method,
            features=features,
            subsample_size=subsample_size,
            with_replacement=sample_with_replacement,
            n_iters=n_iters,
            pairs_only=pairs_only,
            locus=locus,
            seed=seed,
        )
        similarities.append(sim)
    # single comparison, single iteration
    if len(similarities) == 1 and n_iters == 1:
        return similarities[0].similarities[0]
    # single comparison, multiple iterations
    elif len(similarities) == 1:
        return similarities[0]
    # multiple comparisons
    else:
        return RepertoireSimilarities(similarities)


def _calculate_similarity(
    repertoire1: Union[Iterable[Sequence], Iterable[Pair]],
    repertoire2: Union[Iterable[Sequence], Iterable[Pair]],
    name1: str,
    name2: str,
    method: str,
    features: Iterable[str],
    subsample_size: int,
    with_replacement: bool,
    n_iters: int,
    pairs_only: bool,
    locus: Iterable[str],
    seed: Union[int, float, str],
) -> Union[float, RepertoireSimilarity]:
    """
    Calculates the similarity between two repertoires.

    Parameters
    ----------
    repertoire1 : list
        A list of ``abutils.Sequence`` or ``abutils.Pair`` objects.

    repertoire2 : list
        A list of ``abutils.Sequence`` or ``abutils.Pair`` objects.

    name1 : str
        The name of repertoire1.

    name2 : str
        The name of repertoire2.

    method : str
        The similarity method to use.

    features : list
        The features to use for similarity calculation.

    subsample_size : int
        The number of sequences to subsample from each repertoire.

    with_replacement : bool
        Whether to subsample with replacement.

    n_iters : int
        The number of iterations to perform.

    pairs_only : bool
        Whether to use only paired sequences.

    locus : list
        The locus to use for similarity calculation.

    seed : int, float, str
        The seed to use for random number generation.

    Returns
    -------
    RepertoireSimilarity
        A RepertoireSimilarity object.

    """
    # get features
    features1 = _get_features(
        repertoire=repertoire1, features=features, locus=locus, pairs_only=pairs_only
    )
    features2 = _get_features(
        repertoire=repertoire2, features=features, locus=locus, pairs_only=pairs_only
    )
    # seed
    if seed is not None:
        np.random.seed(seed)
    # get similarity function
    similarity_funcs = {
        "morisita-horn": morisita_horn,
        "kullback-leibler": kullback_leibler,
        "jensen-shannon": jensen_shannon,
        "jaccard": jaccard_similarity,
        "bray-curtis": bray_curtis,
        "renkonen": renkonen,
        "cosine": cosine_similarity,
    }
    similarity_func = similarity_funcs.get(method.lower(), None)
    if similarity_func is None:
        err = f"\nThe provided similarity method ({method}) is not supported.\n"
        err += "Supported methods are:\n - "
        err += "\n - ".join(similarity_funcs.keys())
        raise RuntimeError(err)
    # calculate similarity
    similarities = []
    for i in range(n_iters):
        # subsample
        subsample1 = np.random.choice(
            features1, size=subsample_size, replace=with_replacement
        )
        subsample2 = np.random.choice(
            features2, size=subsample_size, replace=with_replacement
        )
        # count features
        counts1 = pd.Series(subsample1).value_counts()
        counts2 = pd.Series(subsample2).value_counts()
        df = pd.DataFrame({"sample1": counts1, "sample2": counts2}).fillna(0)
        # similarity
        similarity = similarity_func(df["sample1"], df["sample2"])
        similarities.append(similarity)
    return RepertoireSimilarity(similarities, method, name1, name2)


def _get_features(
    repertoire: Union[Iterable[Sequence], Iterable[Pair]],
    features: Iterable[str],
    locus: Iterable[str],
    pairs_only: bool = False,
) -> Iterable[str]:
    """
    Gets the requested features from a repertoire.

    Parameters
    ----------
    repertoire : list
        A list of ``abutils.Sequence`` or ``abutils.Pair`` objects.

    features : list
        The features to use for similarity calculation.

    locus : list
        The locus to use for similarity calculation.

    pairs_only : bool
        Whether to use only paired sequences.

    Returns
    -------
    list
        A list of features.

    """
    all_features = []
    # process Pair objects
    if all([isinstance(r, Pair) for r in repertoire]):
        if pairs_only:
            repertoire = [r for r in repertoire if r.is_pair]
        for r in repertoire:
            feat_list = []
            if "IGH" in locus and r.heavy is not None:
                for f in features:
                    feat = r.heavy[f]
                    feat_list.append(str(feat) if feat is not None else "NA")
            if r.light is not None and r.light["locus"] in locus:
                for f in features:
                    feat = r.light[f]
                    feat_list.append(str(feat) if feat is not None else "NA")
            all_features.append("__".join(feat_list))
    else:
        for r in repertoire:
            feat_list = []
            if r["locus"] in locus:
                for f in features:
                    feat = r[f]
                    feat_list.append(str(feat) if feat is not None else "NA")
            all_features.append("__".join(feat_list))
    return all_features


def normalize(
    s1: Iterable[Union[int, float]], s2: Iterable[Union[int, float]]
) -> tuple:
    """
    Normalizes two distributions.

    Parameters
    ----------
    s1 : list
        An iterable of feature frequency values from sample 1.

    s2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    tuple
        A tuple of normalized distributions.

    """
    return s1 / np.sum(s1), s2 / np.sum(s2)


def make_continuous(
    s1: Iterable[Union[int, float]], s2: Iterable[Union[int, float]]
) -> tuple:
    """
    Makes two distributions continuous by removing zero values.

    Parameters
    ----------
    s1 : list
        An iterable of feature frequency values from sample 1.

    s2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    tuple
        A tuple of continuous distributions.

    """
    df = pd.DataFrame({"s1": s1, "s2": s2})
    df = df.fillna(0)
    df = df[(df["s1"] != 0) & (df["s2"] != 0)]
    return df["s1"], df["s2"]


# -----------------------------------
#       SIMILARITY METHODS
# -----------------------------------


def morisita_horn(sample1: Iterable, sample2: Iterable) -> float:
    """
    Calculates the Marista-Horn similarity for two distributions.

    Parameters
    ----------
    sample1 : list
        An iterable of feature frequency values from sample 1.

    sample2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    float
        The Morisita-Horn similarity between sample1 and sample2.

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


def kullback_leibler(sample1: Iterable, sample2: Iterable) -> float:
    """
    Calculates the Kullback-Leibler similarity for two distributions.

    .. note::
        This is not a true metric, but is included for completeness.
        Kullback-Leibler similarity is not symmetric, so the order of the samples
        matters. KL(A, B) != KL(B, A).

    Parameters
    ----------
    sample1 : list
        An iterable of feature frequency values from sample 1.

    sample2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    float
        The Morisita-Horn similarity between sample1 and sample2.

    """
    cont1, cont2 = make_continuous(sample1, sample2)
    norm1, norm2 = normalize(cont1, cont2)
    kl = 0
    for Pi, Qi in zip(norm1, norm2):
        kl += Pi * math.log(Pi / Qi)
    return 1 - kl


def jensen_shannon(sample1: Iterable, sample2: Iterable) -> float:
    """
    Calculates the Jensen-Shannon similarity for two distributions.

    Parameters
    ----------
    sample1 : list
        An iterable of feature frequency values from sample 1.

    sample2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    float
        The Jensen-Shannon similarity between sample1 and sample2.

    """
    norm1, norm2 = normalize(sample1, sample2)
    M = [(n1 + n2) / 2 for n1, n2 in zip(norm1, norm2)]
    js = 0.5 * (kullback_leibler(norm1, M) + kullback_leibler(norm2, M))
    return 1 - js


def jaccard_similarity(sample1: Iterable, sample2: Iterable) -> float:
    """
    Calculates the Jaccard similarity for two distributions.

    Parameters
    ----------
    sample1 : list
        An iterable of feature frequency values from sample 1.

    sample2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    float
        The Jaccard similarity between sample1 and sample2.

    """
    intersection = 0
    union = 0
    for s1, s2 in zip(sample1, sample2):
        intersection += min(s1, s2)
        union += max(s1, s2)
    return intersection / union


def renkonen(sample1: Iterable, sample2: Iterable) -> float:
    """
    Calculates the Renkonen similarity for two distributions.

    Parameters
    ----------
    sample1 : list
        An iterable of feature frequency values from sample 1.

    sample2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    float
        The Renkonen similarity between sample1 and sample2.

    """
    norm1, norm2 = normalize(sample1, sample2)
    renkonen = sum([min(vals) for vals in zip(norm1, norm2)])
    return renkonen


def bray_curtis(sample1: Iterable, sample2: Iterable) -> float:
    """
    Calculates the Bray-Curtis similarity for two distributions.

    Parameters
    ----------
    sample1 : list
        An iterable of feature frequency values from sample 1.

    sample2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    float
        The Bray-Curtis similarity between sample1 and sample2.

    """
    norm1, norm2 = normalize(sample1, sample2)
    sum_diff = sum(abs(x - y) for x, y in zip(norm1, norm2))
    sum_sum = sum(x + y for x, y in zip(norm1, norm2))
    bray_curtis = sum_diff / sum_sum
    return bray_curtis


def cosine_similarity(sample1: Iterable, sample2: Iterable) -> float:
    """
    Calculates the cosine similarity for two distributions.

    Parameters
    ----------
    sample1 : list
        An iterable of feature frequency values from sample 1.

    sample2 : list
        An iterable of feature frequency values from sample 2.

    Returns
    -------
    float
        The cosine similarity between sample1 and sample2.

    """
    numerator = sum(x * y for x, y in zip(sample1, sample2))
    denominator = math.sqrt(sum(x * x for x in sample1)) * math.sqrt(
        sum(y * y for y in sample2)
    )
    cosine = numerator / denominator
    return cosine
