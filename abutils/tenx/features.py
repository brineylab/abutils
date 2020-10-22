#!/usr/bin/env python
# filename: features.py


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
#


__all__ = ['assign_cellhash_groups', 'classify_features',
           'positive_feature_cutoff', 'negative_feature_cutoff']


import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from sklearn.neighbors import KernelDensity

from scipy.signal import argrelextrema
from scipy.stats import scoreatpercentile


def assign_cellhash_groups(df, hashnames=None, rename_dict=None, log_normalize=True,
                           threshold_minimum=4.0, threshold_maximum=10.0, kde_maximum=15.0, 
                           assignments_only=False, return_intermediates=False, debug=False):
    hash_df = df.copy()
    # make sure the dataframe is counts + 1
    if 0 in hash_df.values:
        hash_df += 1
    if log_normalize:
        hash_df = hash_df.applymap(np.log2)
    if hashnames is None:
        hashnames = hash_df.columns.values
    if rename_dict is None:
        rename_dict = {h: h for h in hashnames}

    thresholds = {}
    for hashname in hashnames:
        if debug:
            print(hashname)
        thresholds[hashname] = positive_feature_cutoff(hash_df[hashname],
                                                       threshold_minimum=threshold_minimum,
                                                       threshold_maximum=threshold_maximum,
                                                       kde_maximum=kde_maximum,
                                                       debug=debug)
    
    if debug:
        print('THRESHOLDS')
        print('----------')
        for hashname in hashnames:
            print(f'{hashname}: {thresholds[hashname]}')

    assignments = []
    for i, row in hash_df[hashnames].iterrows():
        a = [h for h in hashnames if row[h] >= thresholds[h]]
        if len(a) == 1:
            assignment = rename_dict.get(a[0], 'not found')
        elif len(a) > 1:
            assignment = 'doublet'
        else:
            assignment = 'unassigned'
        assignments.append(assignment)
    hash_df['assignment'] = assignments
    
    if assignments_only:
        return hash_df['assignment']
    else:
        if return_intermediates:
            return hash_df
        else:
            return hash_df[hashnames + ['assignment']]


def classify_features(df, positive=None, negative=None,
                      positive_threshold_maximum=10.0, positive_threshold_minimum=4.0, positive_kde_maximum=15.0,
                      negative_threshold_maximum=10.0, negative_threshold_minimum=4.0, negative_kde_maximum=15.0,
                      negative_threshold_denominator=2, classifications_only=False, debug=False):
    '''
    Docstring for classify_features.
    '''
    df = df.copy()
    data = {}
    
    if positive is None:
        positive = []
    if negative is None:
        negative = []
    
    for p in positive:
        if debug:
            print(p)
        cutoff = positive_feature_cutoff(df[p],
                                         threshold_maximum=positive_threshold_maximum,
                                         threshold_minimum=positive_threshold_minimum,
                                         kde_maximum = positive_kde_maximum,
                                         debug=debug)
        data[p] = df[p] >= cutoff
    for n in negative:
        if debug:
            print(n)
        cutoff = negative_feature_cutoff(df[n],
                                         threshold_maximum=negative_threshold_maximum,
                                         threshold_minimum=negative_threshold_minimum,
                                         kde_maximum = negative_kde_maximum,
                                         denominator=negative_threshold_denominator,
                                         debug=debug)
        data[n] = df[n] > cutoff
    
    if classifications_only:
        return pd.DataFrame(data)
    else:
        data = {'{}_classified'.format(k): v for k, v in data.items()}
        return df.join(pd.DataFrame(data))
        
    


def positive_feature_cutoff(vals, threshold_maximum=10.0, threshold_minimum=4.0,
                            kde_maximum=15.0, debug=False):
    a = np.array(vals)
    k = _bw_silverman(a)
    kde = KernelDensity(kernel='gaussian', bandwidth=k).fit(a.reshape(-1, 1))
    s = np.linspace(0, kde_maximum, num=int(kde_maximum * 100))
    e = kde.score_samples(s.reshape(-1,1))
    
    all_min, all_max = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
    if len(all_min) > 1:
        _all_min = np.array([m for m in all_min if s[m] <= threshold_maximum and s[m] >= threshold_minimum])
        min_vals = zip(_all_min, e[_all_min])
        mi = sorted(min_vals, key=lambda x: x[1])[0][0]
        cutoff = s[mi]
    elif len(all_min) == 1:
        mi = all_min[0]
        cutoff = s[mi]
    else:
        cutoff = None

    if debug:
        plt.plot(s, e)
        plt.show()
        print('bandwidth: {}'.format(k))
        print('local minima: {}'.format(s[all_min]))
        print('local maxima: {}'.format(s[all_max]))
        if cutoff is not None:
            print('cutoff: {}'.format(cutoff))
        else:
            print('WARNING: no local minima were found, so the threshold could not be calculated.')
        print('\n\n')

    return cutoff


def negative_feature_cutoff(vals, threshold_maximum=10.0, threshold_minimum=4.0,
                            kde_maximum=15.0, denominator=2, debug=False):
    a = np.array(vals)
    k = _bw_silverman(a)
    kde = KernelDensity(kernel='gaussian', bandwidth=k).fit(a.reshape(-1, 1))
    s = np.linspace(0, kde_maximum, num=int(kde_maximum * 100))
    e = kde.score_samples(s.reshape(-1,1))
    
    all_min, all_max = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]
    if len(all_min) > 1:
        _all_min = np.array([m for m in all_min if s[m] <= threshold_maximum and s[m] >= threshold_minimum])
        min_vals = zip(_all_min, e[_all_min])
        mi = sorted(min_vals, key=lambda x: x[1])[0][0]
        i = list(all_min).index(mi)
        try:
            if s[mi] > s[all_max[i]]:
                ma = all_max[i]
            else:
                ma = all_max[i - 1]
        except IndexError:
            ma = all_max[i-1]
        cutoff = (s[mi] + s[ma]) / denominator
    elif len(all_min) == 1:
        mi = all_min[0]
        ma = all_max[0]
        cutoff = (s[mi] + s[ma]) / denominator
    else:
        cutoff = None

    if debug:
        plt.plot(s, e)
        plt.show()
        print('bandwidth: {}'.format(k))
        print('local minima: {}'.format(s[all_min]))
        print('local maxima: {}'.format(s[all_max]))
        if cutoff is not None:
            print('cutoff: {}'.format(cutoff))
        else:
            print('WARNING: no local minima were found, so the threshold could not be calculated.')
        print('\n\n')

    return cutoff
    
    
def _bw_silverman(x):
    normalize = 1.349
    IQR = (scoreatpercentile(x, 75) - scoreatpercentile(x, 25)) / normalize
    std_dev = np.std(x, axis=0, ddof=1)
    if IQR > 0:
        A = np.minimum(std_dev, IQR)
    else:
        A = std_dev
    n = len(x)
    return .9 * A * n ** (-0.2)




