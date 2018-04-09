#!/usr/bin/env python
# filename: summary.py


#
# Copyright (c) 2017 Bryan Briney
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


from __future__ import print_function

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from abstar.core.germline import get_germlines


__all__ = ['germline_use_plot', 'cdr3_length_plot', 'vj_heatmap']


#==============================
#           PLOTS
#==============================


def germline_use_plot(seqs, fig_file=None, level='gene', gene='V', species='human', chain='heavy'):
    from .base import barplot
    seqs = [s for s in seqs if s['chain'] == chain]
    germs = get_germlines(species, gene, chain=chain)
    if level == 'fam':
        keys = list(set([g.name.split('-')[0] for g in germs]))
        size = (6, 4)
    elif level == 'gene':
        keys = list(set([g.name.split('*')[0] for g in germs]))
        if gene == 'J':
            size = (6, 4)
        elif gene == 'D':
            size = (8, 4)
        else:
            size = (10, 4)
    gtype = '{}_gene'.format(gene.lower())
    short_chain = chain[0].upper()
    data = [s[gtype][level].rstrip('D') for s in seqs if gtype in s]
    x, y = _aggregate(data, keys=keys)
    x = [i.replace('IG{}{}'.format(short_chain, gene), '{}{}'.format(gene, short_chain)) for i in x]
    colors = _get_germline_plot_colors(x)
    barplot(x, y,
            colors,
            fig_file=fig_file,
            ylabel='Frequency (%)',
            rotate_xtick_labels=True,
            size=size)


def cdr3_length_plot(seqs, fig_file=None, max_len=40, chain='heavy', color=None):
    from .base import barplot
    if chain == 'light':
        chain = ['kappa', 'lambda']
    else:
        chain = [chain, ]
    seqs = [s for s in seqs if s['chain'] in chain]
    cdr3s = [s['cdr3_len'] for s in seqs if s['cdr3_len'] > 0 and s['cdr3_len'] <= max_len]
    x, y = _aggregate(cdr3s, keys=list(range(1, max_len + 1)))
    color = color if color is not None else sns.hls_palette(7)[4]
    x_title = 'CDR3 Length (AA)'
    y_title = 'Frequency (%)'
    size = ((max_len / 5.0) + 1, 4)
    barplot([str(i) for i in x], y,
            color,
            fig_file=fig_file,
            xlabel=x_title,
            ylabel=y_title,
            grid=True,
            size=size,
            xfontsize=7)


def vj_heatmap(seqs, fig_file=None, species='human', chain='heavy'):
    from .base import heatmap
    seqs = [s for s in seqs if s['chain'] == chain]
    vj_data = _group_by_vj(seqs, species, chain)
    vj_df = pd.DataFrame(vj_data)
    heatmap(vj_df.transpose(), fig_file=fig_file)



#==============================
#        PLOT UTILS
#==============================


def _aggregate(data, norm=True, sort_by='value', keys=None):
    '''
    Counts the number of occurances of each item in 'data'.

    Inputs
    data: a list of values.
    norm: normalize the resulting counts (as percent)
    sort_by: how to sort the retured data. Options are 'value' and 'count'.

    Output
    a non-redundant list of values (from 'data') and a list of counts.
    '''
    if keys:
        vdict = {k: 0 for k in keys}
        for d in data:
            if d in keys:
                vdict[d] += 1
    else:
        vdict = {}
        for d in data:
            vdict[d] = vdict[d] + 1 if d in vdict else 1
    vals = [(k, v) for k, v in vdict.items()]
    if sort_by == 'value':
        vals.sort(key=lambda x: x[0])
    else:
        vals.sort(key=lambda x: x[1])
    xs = [v[0] for v in vals]
    if norm:
        raw_y = [v[1] for v in vals]
        total_y = sum(raw_y)
        ys = [100. * y / total_y for y in raw_y]
    else:
        ys = [v[1] for v in vals]
    return xs, ys


def _get_germline_plot_colors(data):
    fams = [d.split('-')[0].replace('/OR', '') for d in data]
    nr_fams = sorted(list(set(fams)))
    num_colors = len(nr_fams)
    rgbs = sns.hls_palette(num_colors)
    rgb_dict = {i[0]: i[1] for i in zip(nr_fams, rgbs)}
    return [rgb_dict[f] for f in fams]


def _group_by_vj(data, species, chain):
    # from abtools.utils.germlines import germlines
    vs = list(set([g.name.split('*')[0] for g in get_germlines(species, 'V', chain=chain)]))
    js = list(set([g.name.split('*')[0] for g in get_germlines(species, 'J', chain=chain)]))
    vj = {}
    for v in vs:
        vj[v] = {}
        for j in js:
            vj[v][j] = 0
    for d in data:
        v = d['v_gene']['gene']
        j = d['j_gene']['gene'].rstrip('P')
        if v[:3] != j[:3]:
            continue
        if v not in vj:
            continue
        if j not in vj[v]:
            continue
        vj[v][j] = vj[v][j] + 1 if j in vj[v] else 1
    total = len(data)
    for v in list(vj.keys()):
        for j in list(vj[v].keys()):
            vj[v][j] = 100. * vj[v][j] / total
    return vj

