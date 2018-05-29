#!/usr/bin/env python
# filename: plots.py


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


__all__ = ['barplot', 'heatmap']


def barplot(x, y, colors, fig_file=None, xlabel=None, ylabel=None, rotate_xtick_labels=False, grid=False, size=None, xfontsize=None):
    sns.set_style('ticks')
    # set bar locations and width
    ind = np.arange(len(x))
    width = 0.75
    # plot objects
    if size:
        fig = plt.figure(figsize=size)
    else:
        fig = plt.figure()
    ax = fig.add_subplot(111)
    # axis limits and ticks
    ax.set_ylim(0, 1.05 * max(y))
    ax.set_xlim(-width, len(ind) - (width / 2))
    # ax.set_xticks(ind + width / 2)
    ax.set_xticks(ind)
    xtick_names = ax.set_xticklabels(x)
    if rotate_xtick_labels:
        plt.setp(xtick_names, rotation=90, fontsize=7)
    if grid:
        ax.yaxis.grid(True, alpha=0.5)
    # axis labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xfontsize is not None:
        plt.setp(xtick_names, fontsize=xfontsize)
    ax.tick_params(axis='x', which='both', top='off', length=3, pad=1.5)
    ax.tick_params(axis='y', which='both', right='off', length=3, pad=1.5)
    # make the plot
    bar = ax.bar(ind, y, width, color=colors)
    fig.tight_layout()
    if fig_file is not None:
        plt.savefig(fig_file)
    else:
        plt.show()
    plt.close()



def heatmap(df, fig_file=None):
    sns.set()
    # set up plot, determine plot size
    h, w = df.shape
    f, ax = plt.subplots(figsize=(w / 1.75, h / 3))
    sns.heatmap(df,
                square=True,
                cmap='Blues',
                cbar=True,
                cbar_kws={'orientation': 'horizontal',
                          'fraction': 0.02,
                          'pad': 0.02,
                          'shrink': 0.675})
    # adjust labels
    ax.xaxis.tick_top()
    plt.xticks(rotation=90)
    f.tight_layout()
    # make the plot
    if fig_file is not None:
        plt.savefig(fig_file)
    else:
        plt.show()
    plt.close()