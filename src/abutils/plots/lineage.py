#!/usr/bin/env python
# filename: lineage.py


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


import colorsys
import random

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from natsort import natsorted

from ..utils.color import get_cmap


__all__ = ['lineage_donut', 'annotated_donut']


def lineage_donut(lineages, figfile=None, figsize=(6, 6), pairs_only=False,
                  monochrome_color=None, singleton_color='lightgrey', shuffle_colors=False, seed=1234,
                  text_kws={}, pie_kws={}, fontsize=28, linewidth=2):
    lineages = sorted(lineages, key=lambda x: x.size(pairs_only), reverse=True)
    non_singletons = [l for l in lineages if l.size(pairs_only) > 1]
    singleton_count = sum([1 for l in lineages if l.size(pairs_only) == 1])
    lineage_sizes = [l.size(pairs_only) for l in lineages if l.size(pairs_only) > 1] + [singleton_count]
    if monochrome_color is not None:
        colors = _get_monochrome_colors(monochrome_color, len(non_singletons))
        # we shuffle the colors differently if we're using a monochrome palette, because
        # we want to have the first color (largest lineage) always be the user-supplied
        # monochrome_color. We only want to shuffle the colors starting with the second one.
        if shuffle_colors:
            primary = colors[0]
            secondary = colors[1:]
            random.seed(seed)
            random.shuffle(secondary)
            colors = [primary] + secondary
    else:
        colors = _get_donut_colors(len(non_singletons))
        if shuffle_colors:
            random.seed(seed)
            random.shuffle(colors)
    colors += [singleton_color]

    plt.figure(figsize=figsize)
    ax = plt.gca()
    # fig, ax = plt.subplots()
    ax.axis('equal')
    width = 0.55
    kwargs = dict(colors=colors, startangle=90)
    for k, v in pie_kws.items():
        kwargs[k] = v
    inside, _ = ax.pie(lineage_sizes, radius=1, pctdistance=1 - width / 2, **kwargs)
    plt.setp(inside, width=width, edgecolor='white')

    for w in inside:
        w.set_linewidth(linewidth)

    kwargs = dict(size=fontsize, color='k', va='center', fontweight='bold')
    for k, v in text_kws.items():
        kwargs[k] = v
    ax.text(0, 0, str(sum(lineage_sizes)), ha='center', **kwargs)

    plt.tight_layout()

    if figfile is not None:
        plt.savefig(figfile)
    else:
        plt.show()


def annotated_donut(lineages, annotations, annotation_colors=None, annotation_order=None, annotation_fontsize=20,
                    unannotated_color='grey', unannotated_label='unknown', arc_width=0.125, donut_width=0.425,
                    count_fontsize=48, linewidth=2, arc_linewidth=15, figfile=None, figsize=(8, 8), pairs_only=False,
                    monochrome_color=None, singleton_color='lightgrey', shuffle_colors=False, seed=1234,
                    text_kws={}, pie_kws={}):

    all_sizes = []
    all_colors = []
    if annotation_order is None:
        annotation_order = list(natsorted(annotations.keys()))
    if annotation_colors is None:
        annotation_colors = {a: c for a, c in zip(annotation_order, _get_donut_colors(len(annotation_order)))}

    # annotated_lineages
    for ann in annotation_order:
        alineages = [l for l in lineages if l.name in annotations[ann]]
        non_singletons = [l for l in alineages if l.size(pairs_only) > 1]
        singleton_count = sum([1 for l in alineages if l.size(pairs_only) == 1])
        lineage_sizes = [l.size(pairs_only) for l in alineages if l.size(pairs_only) > 1] + [singleton_count]
        colors = _get_monochrome_colors(annotation_colors[ann], len(non_singletons))
        if shuffle_colors:
            primary = colors[0]
            secondary = colors[1:]
            random.seed(seed)
            random.shuffle(secondary)
            colors = [primary] + secondary
        colors += [singleton_color]
        all_sizes.append(np.array(lineage_sizes))
        all_colors.append(colors)

    # unannotated_lineages
    ulineages = [l for l in lineages if all([l.name not in names for names in annotations.values()])]
    non_singletons = [l for l in ulineages if l.size(pairs_only) > 1]
    singleton_count = sum([1 for l in ulineages if l.size(pairs_only) == 1])
    lineage_sizes = [l.size(pairs_only) for l in ulineages if l.size(pairs_only) > 1] + [singleton_count]
    colors = _get_monochrome_colors(unannotated_color, len(non_singletons))
    if shuffle_colors:
        primary = colors[0]
        secondary = colors[1:]
        random.seed(seed)
        random.shuffle(secondary)
        colors = [primary] + secondary
    colors += [singleton_color]
    all_sizes.append(np.array(lineage_sizes))
    all_colors.append(colors)
    
    # sizes and colors
    sizes = np.array(all_sizes)
    outer_sizes = [s.sum() for s in sizes]
    inner_sizes = [s for sublist in sizes for s in sublist]
    colors = np.array(all_colors)
    outer_colors = [annotation_colors[a] for a in annotation_order] + [unannotated_color]
    inner_colors = [c for sublist in colors for c in sublist]
    
    # make the figure
    plt.subplots(figsize=figsize, subplot_kw=dict(aspect="equal"))
    ax = plt.gca()
    # annotation arc
    wedges, texts = ax.pie(outer_sizes, radius=1, colors=outer_colors, startangle=90,
                           wedgeprops=dict(width=arc_width, edgecolor='w', linewidth=arc_linewidth))
    # donut plot
    ax.pie(inner_sizes, radius=1-arc_width, colors=inner_colors, startangle=90,
           wedgeprops=dict(width=donut_width, edgecolor='w'))
    # sequence count
    kwargs = dict(size=count_fontsize, color='k', va='center', fontweight='bold')
    for k, v in text_kws.items():
        kwargs[k] = v
    ax.text(0, 0, str(sum(inner_sizes)), ha='center', **kwargs)
    
    # add annotation labels
    annot_labels = annotation_order + [unannotated_label]
    kw = dict(arrowprops=dict(arrowstyle="-"),
              fontsize=annotation_fontsize,
              fontweight='medium',
              zorder=0, va="center")
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        kw['arrowprops'].update({'color': outer_colors[i]})
        ax.annotate(annot_labels[i], xy=(x, y),
                    xytext=(1.15*np.sign(x), 1.2*y),
                    color=outer_colors[i],
                    horizontalalignment=horizontalalignment, **kw)
    
    # show or save
    plt.tight_layout()
    if figfile is not None:
        plt.savefig(figfile)
    else:
        plt.show()








def _get_donut_colors(N):
    HSV_tuples = [(x * 1.0 / N, 0.8, 0.9) for x in range(N)]
    RGB_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))[::-1]
    return RGB_tuples


def _get_monochrome_colors(monochrome_color, N):
    cmap = get_cmap(from_color=monochrome_color)
    # this is a bit convoluted, but what's happening is we're getting different colormap
    # values (which range from 0 to 1). Calling cmap(i) returns an rgba tuple, but we just need
    # the rbg, so we drop the a. To make sure that one of the colors isn't pure white,
    # we ask np.linspace() for one more value than we need, reverse the list of RGB tuples
    # so that it goes from dark to light, and drop the lightest value
    RGB_tuples = [cmap(i)[:-1] for i in np.linspace(0, 1, N + 1)][::-1][:-1]
    return RGB_tuples