#!/usr/bin/env python
# filename: pl.py


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


import sys

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import seaborn as sns

from natsort import natsorted



def feature_kde(data, x, y, hue=None, hue_order=None, colors=None, thresh=0.1,
                show_scatter=True, scatter_size=5, scatter_alpha=0.2,
                fill=False, kde_fill_alpha=0.7, kde_line_alpha=1.0,
                highlight_index=None, highlight_x=None, highlight_y=None, highlight_marker='x',
                highlight_size=90, highlight_color='k', highlight_name=None, highlight_alpha=0.8,
                xlabel=None, ylabel=None, equal_axes=True,
                legend_loc='best', legend_title=None, show_legend_title=True,
                return_ax=False, figsize=[6, 6], figfile=None):
    '''
    Produces a 2-dimensional KDE plot of two features.

    Args:

        data (pd.DataFramne): A pandas dataframe containing the input data.

        x (str): Name of the column in ``df`` containing the feature to be plotted on the x-axis.

        y (str): Name of the column in ``df`` containing the feature to be plotted on the y-axis.

        hue (str): Name of the column in ``df`` containing categories for hue values. For scatter plots, 
                   the categories in ``hue`` will be plotted as different colored points. For KDE plots,
                   ``hue```` categories will each be plotted as differently colored KDE plots
                   on the same plot. 

        hue_order (iterable): Iterable of hue categories, in the order they should be plotted and listed
                              in the legend. If ```hue_order``` contains only a subset of the categories
                              present in ```data[hue]```, only the categories supplied in ```hue_order```
                              will be plotted.

        colors (iterable): List of colors to be used for ```'hue'``` categories. If ```'colors'``` is
                           shorter than the list of hue categories, colors will be reused.

        thresh (float): Threshold for the KDE. Default is ```0.1```.
        
        show_scatter (bool): Show the scatterplot beneath the transparent KDE plot. Default is ```True```.

        scatter_size (int, float): Size of the scatter points. Default is ```5```.

        scatter_alpha (float): Alpha of the scatter points. Default is ```0.2```.

        fill (bool): Fill the KDE plot. Default is ```True```.

        kde_fill_alpha (float): Alpha for the filled KDE plot. If ```fill``` is ```False```,
                                this option is ignored. Default is ```0.7```.
        
        kde_line_alpha (float): Alpha for the KDE plot lines. Default is ```1.0```.

        highlight_index (iterable): An iterabile of index names (present in ```data```) of points
                                    to be highlighted on the KDE plot. If provided, ```highlight_x```
                                    and ```highlight_y``` are ignored.

        highlight_x (iterable): An iterable of x-values for highlighted points. Also requires
                                ```highlight_y```.
        
        highlight_y (iterable): An iterable of y-values for highlighted points. Also requires
                                ```highlight_x```.
        
        highlight_marker (str): The marker style to be used for highlight points. Accepts 
                                standard matplotlib marker styles. Default is ```'x'```. 

        highlight_size (int): Size of the highlight marker. Default is ```90```.

        highlight_color (string or RGB list): Color of the highlight points. Default is black.

        highlight_name (str): Name of the highlights, to be used in the legend. If not supplied,
                              highlight points will not be included in the legend.
        
        highlight_alpha (float): Alpha for the highlight points. Default is ```0.8```.

        xlabel (str): Label for the x-axis. By default, the value for ```x``` is used.

        ylabel (str): Label for the y-axis. By default, the value for ```y``` is used.

        equal_axes (bool): If ```True```, the the limits of the x- and y-axis will be equal.
                           Default is ```True```.
        
        legend_loc (str): Location for the legend. Uses standard matplotlib locations. Default
                          is ```'best'```.

        legend title (str): Title for the legend. By default, ```hue``` is used.

        show_legend_title (bool): Whether or not to show the legend title. Default is ```True```.

        return_ax (bool): If ```True```, return the plot's ```ax``` object. Will not show or save
                          the plot. Default is ```False```.

        figsize (list): A list containg the dimensions of the plot. Default is ```[6, 6]```.

        figfile (str): Path to which the figure will be saved. If not provided, the figure will be
                       shown but not saved to file.
    '''

    # input data
    data = data.copy()
    if x not in data.columns.values:
        print('"{}" is not a column in the supplied dataframe'.format(x))
        return
    if y not in data.columns.values:
        print('"{}" is not a column in the supplied dataframe'.format(y))
        return

    # hue
    if hue is not None:
        if hue_order is None:
            hue_order = natsorted(list(set(data[hue])))
        data = data[data[hue].isin(hue_order)]
    else:
        hue_order = []

    # colors
    n_colors = max(1, len(hue_order))
    if colors is None:
        colors = sns.color_palette(n_colors=n_colors)
        
    plt.figure(figsize=figsize)

    # scatterplots
    if show_scatter:
        if hue_order:
            for h, c in zip(hue_order, colors):
                d = data[data[hue] == h]
                plt.scatter(d[x], d[y], c=[c], s=scatter_size,
                            alpha=scatter_alpha, linewidths=0)
        else:
            plt.scatter(data[x], data[y], c=[colors[0]], s=scatter_size,
                            alpha=scatter_alpha, linewidths=0)

    # kdeplot
    if fill:
        sns.kdeplot(data=data, x=x, y=y, hue=hue, fill=True, alpha=kde_fill_alpha,
                    hue_order=hue_order, palette=colors, thresh=thresh)
    ax = sns.kdeplot(data=data, x=x, y=y, hue=hue, alpha=kde_line_alpha,
                     hue_order=hue_order, palette=colors, thresh=thresh)
    
    # highlighted points
    highlight = any([highlight_index is not None, all([highlight_x is not None, highlight_y is not None])])
    if highlight:
        if highlight_index is not None:
            hidata = data.loc[highlight_index]
            highlight_x = hidata[x]
            highlight_y = hidata[y]
        plt.scatter(highlight_x, highlight_y, zorder=10,
                    s=highlight_size,
                    c=highlight_color,
                    alpha=highlight_alpha,
                    marker=highlight_marker)

    # legend
    legend_labels = hue_order
    if fill:
        handles = []
        for c in colors:
            f = Patch(fc=c, alpha=kde_fill_alpha / 3)
            e = Patch(ec=c, fill=False, lw=1.5)
            handles.append((f, e))
#         handles = [Patch(fc=c, ec=c, alpha=kde_fill_alpha / 3, label=h) for c, h in zip(colors, hue_order)]
    else:
        handles = [Line2D([0], [0], color=c, label=h) for c in colors]
    if highlight_name is not None:
        legend_labels.append(highlight_name)
        handles.append(Line2D([0], [0], marker=highlight_marker, color='w',
                              mec=highlight_color,
                              mfc=highlight_color,
                              ms=highlight_size / 10))
    if show_legend_title:
        legend_title = legend_title if legend_title is not None else hue
    else:
        legend_title = None
    ax.legend(handles, legend_labels, loc=legend_loc, fontsize=13, title=legend_title, frameon=False)
    
    # style the plot
    ax.set_xlabel(xlabel if xlabel is not None else x, fontsize=16)
    ax.set_ylabel(ylabel if ylabel is not None else y, fontsize=16)
    ax.tick_params(axis='both', labelsize=13)

    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_position(('outward', 10))

    if equal_axes:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        axlim = [min([xlim[0], ylim[0]]), max([xlim[1], ylim[1]])]
        ax.set_xlim(axlim)
        ax.set_ylim(axlim)
    
    if return_ax:
        return ax
    elif figfile is not None:
        plt.tight_layout()
        plt.savefig(figfile)
    else:
        plt.show()


    
def feature_scatter(data, x, y, hue=None, hue_order=None, colors=None, marker='o', size=20, alpha=0.6,
                    highlight_index=None, highlight_x=None, highlight_y=None, highlight_marker='o',
                    highlight_size=90, highlight_color='k', highlight_name=None, highlight_alpha=0.9,
                    xlabel=None, ylabel=None, equal_axes=True,
                    legend_loc='best', legend_title=None, show_legend_title=True,
                    return_ax=False, figsize=[6, 6], figfile=None):
    '''
    Docstring for feature_scatter.
    '''
    # input data
    data = data.copy()
    if x not in data.columns.values:
        print('"{}" is not a column in the supplied dataframe'.format(x))
        return
    if y not in data.columns.values:
        print('"{}" is not a column in the supplied dataframe'.format(y))
        return
    
    # hue
    if hue is not None:
        if hue_order is None:
            hue_order = natsorted(list(set(data[hue])))
        data = data[data[hue].isin(hue_order)]
    else:
        hue_order = []

    # colors
    n_colors = max(1, len(hue_order))
    if colors is None:
        colors = sns.color_palette(n_colors=n_colors)
        
    plt.figure(figsize=figsize)
    ax = plt.gca()

    # scatterplot
    if hue_order:
        for h, c in zip(hue_order, colors):
            d = data[data[hue] == h]
            plt.scatter(d[x], d[y], c=[c], s=size, marker=marker,
                        alpha=alpha, linewidths=0, label=h)
    else:
        plt.scatter(data[x], data[y], c=[colors[0]], s=size, marker=marker,
                        alpha=alpha, linewidths=0, label=h)

    # highlighted points
    highlight = any([highlight_index is not None, all([highlight_x is not None, highlight_y is not None])])
    if highlight:
        if highlight_index is not None:
            hidata = data.loc[highlight_index]
            highlight_x = hidata[x]
            highlight_y = hidata[y]
        plt.scatter(highlight_x, highlight_y, zorder=10,
                    s=highlight_size,
                    c=highlight_color,
                    alpha=highlight_alpha,
                    marker=highlight_marker,
                    label=highlight_name)
    # legend
    ax.legend(loc=legend_loc, fontsize=13, title=legend_title, frameon=False)
    
    # style the plot
    ax.set_xlabel(xlabel if xlabel is not None else x, fontsize=16)
    ax.set_ylabel(ylabel if ylabel is not None else y, fontsize=16)
    ax.tick_params(axis='both', labelsize=13)

    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        ax.spines[spine].set_position(('outward', 10))

    if equal_axes:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        axlim = [min([xlim[0], ylim[0]]), max([xlim[1], ylim[1]])]
        ax.set_xlim(axlim)
        ax.set_ylim(axlim)
    
    if return_ax:
        return ax
    elif figfile is not None:
        plt.tight_layout()
        plt.savefig(figfile)
    else:
        plt.show()

