#!/usr/bin/env python
# filename: bar.py


#
# Copyright (c) 2022 Bryan Briney
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


from collections import Counter
import itertools
import os
import re
import sys

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import scanpy as sc

import seaborn as sns

from anndata import AnnData

from natsort import natsorted

from abutils.utils.color import get_cmap
from abutils.utils.utilities import nested_dict_lookup


def bar(
    x=None,
    y=None,
    hue=None,
    data=None,
    order=None,
    hue_order=None,
    palette=None,
    color=None,
    alt_color="#D3D3D3",
    normalize=False,
    highlight=None,
    highlight_color=None,
    orientation="vertical",
    plot_kwargs=None,
    legend_kwargs=None,
    hide_legend=False,
    xlabel=None,
    ylabel=None,
    xlabel_fontsize=16,
    ylabel_fontsize=16,
    xtick_labelsize=14,
    ytick_labelsize=14,
    xtick_labelrotation=0,
    ytick_labelrotation=0,
    show=False,
    figsize=None,
    figfile=None,
):
    """
    Produces a bar plot of categorical data. For data with distinct batches, a stacked 
    bar plot will be constructed.

    Parameters
    ----------

    x : str or list, optional  
        Name of a column in `data` or an iterable of values to be plotted on the 
        x-axis. Required if `orientation` is ``'horizontal'``.

    y : str or list, optional  
        Name of a column in `data` or an iterable of values to be plotted on the 
        y-axis. Required if `orientation` is ``'horizontal'``.

    hue : str, optional  
        Name of a column in `data` or an iterable of hue categories to be used to 
        group data into stacked bars. If not provided, an un-stacked bar plot is created.  

    data : pandas.DataFrame, optional  
        A ``DataFrame`` object containing the input data. If provided, `x` and/or `y` should 
        be column names in `data`.
        
    order : iterable object, optional  
        List of `x` or `y` categories in the order they should be plotted. If `order` contains a 
        subset of all categories found in `x` or `y`, only the supplied categories will be plotted.  
        If not provided, categories will be plotted in ``natsort.natsorted()`` order.
        
    hue_order : iterable object, optional  
        List of `hue` categories in the order they should be plotted. If `hue_order` contains a 
        subset of all categories found in `hue`, only the supplied categories will be plotted.  
        If not provided, `hue` categories will be plotted in ``natsort.natsorted()`` order.
        
    palette : dict, optional  
        Dictionary mapping `hue`, `x` or `y` names to colors. If if keys in `palette` match 
        more than one category, `hue` categories take priority. If `palette` is not provided, 
        bars are colored using `color` (if `hue` is ``None``) or a palette is generated 
        automatically using ``sns.hls_palette()``.   

    color : str or iterable, optional  
        Single color to be used for the bar plot. If not provided, the first color in the 
        default ``Seaborn`` color palette will be used. If `highlight` is provided but 
        `highlight_color` is not, `color` will be used to color highlighted bars.  
        
    alt_color : str or iterable, default='#D3D3D3'  
        Alternate color for the bar plot. Used to color categories not provided in `palette` 
        or to color categories not present in `highlight`.  

    orientation : str, optional
        Orientation of the plot. Options are ``'vertical'`` or ``'horizontal'``. Default is 
        ``'vertical'``.
        
    normalize : bool, default=False  
        If ``True``, normalized frequencies are plotted instead of raw counts. If multiple `hue`
        categories are present, data will be normalized such that all 
        bars extend from [0,1] and each stacked bar is sized according to the `hue`'s  fraction. 
        If `hue` is not provided or there is only one `hue` category, the entire 
        dataset is normalized.
        
    highlight : iterable, optional
        List of `x` or `hue` categories to be highlighted. If `highlight_color` is provided, 
        categories in `highlight` will use `highlight_color` and all others will use `alt_color`. 
        If `highlight_color` is not provided, `palette` will be used. If both `highlight_color` 
        and `palette` are not provided, `color` will be used. 
        
    highlight_color : str or iterable, optional  
        Color to be used for categories in `highlight`. If  

    plot_kwargs : dict, optional  
        Dictionary containing keyword arguments that will be passed to ``pyplot.bar()``.

    legend_kwargs : dict, optional  
        Dictionary containing keyword arguments that will be passed to ``ax.legend()``.

    hide_legend : bool, default=False  
        By default, a plot legend will be shown if multiple batches are plotted. If ``True``, 
        the legend will not be shown.  
        
    xlabel : str, optional  
        Text for the X-axis label. 

    ylabel : str, optional  
        Text for the Y-axis label.  
        
    xlabel_fontsize : int or float, default=16  
        Fontsize for the X-axis label text.

    ylabel_fontsize : int or float, default=16  
        Fontsize for the Y-axis label text.

    xtick_labelsize : int or float, default=14  
        Fontsize for the X-axis tick labels.  

    ytick_labelsize : int or float, default=14  
        Fontsize for the Y-axis tick labels.  

    xtick_labelrotation : int or float, default=0  
        Rotation of the X-axis tick labels.  
    
    ytick_labelrotation : int or float, default=0  
        Rotation of the Y-axis tick labels. 

    show :bool, default=False  
        If ``True``, plot is shown and the plot ``Axes`` object is not returned. Default
        is ``False``, which does not call ``pyplot.show()`` and returns the ``Axes`` object.

    figsize : iterable object, default=[6, 4]  
        List containing the figure size (as ``[x-dimension, y-dimension]``) in inches. 

    figfile : str, optional  
        Path at which to save the figure file. If not provided, the figure is not saved
        and is either shown (if `show` is ``True``) or the ``Axes`` object is returned.  
    """
    # process input data
    if orientation == "horizontal":
        if y is not None:
            x, y = y, x
    if data is None:
        _data = {}
        if x is not None:
            _data["x"] = x
            x = "x"
        if y is not None:
            _data["y"] = y
            y = "y"
        if hue is not None:
            _data["hue"] = hue
            hue = "hue"
        data = pd.DataFrame(_data)

    # figure size
    if figsize is None:
        if orientation == "horizontal":
            figsize = [4, 6]
        else:
            figsize = [6, 4]

    if order is None:
        x_vals = natsorted(data[x].unique())
    else:
        x_vals = order

    # process hue, if provided
    if hue is not None:
        hue_vals = data[hue]
        hue_order = hue_order if hue_order is not None else natsorted(set(hue_vals))
        hue_batches = [
            data[[_hue_val == h for _hue_val in hue_vals]] for h in hue_order
        ]
    else:
        hue_order = [
            None,
        ]
        hue_batches = [
            data,
        ]

    # process batches
    batch_data = []
    for batch in hue_batches:
        if y is not None:
            y_dict = {}
            for _x, _y in zip(batch[x], batch[y]):
                if _x in y_dict:
                    y_dict[_x] += float(_y)
                else:
                    y_dict[_x] = float(_y)
            batch_data.append(y_dict)
        else:
            y_dict = Counter(batch[x])
            batch_data.append(y_dict)
    if normalize:
        if len(batch_data) > 1:
            ytots = {xval: sum([b[xval] for b in batch_data]) for xval in x_vals}
            for y_dict in batch_data:
                for xval in x_vals:
                    yval = y_dict.get(xval, 0)
                    y_dict[xval] = yval / ytots[xval]
        else:
            for ydict in batch_data:
                tot = sum(ydict.values())
                for xval in x_vals:
                    yval = y_dict.get(xval, 0)
                    y_dict[xval] = yval / tot

    # colors
    if palette is None:
        if len(hue_batches) > 1:
            palette = {h: c for h, c in zip(hue_order, sns.hls_palette(len(hue_order)))}
        else:
            palette = {}
        # palette is a dict assigning colors to x and/or hue categories
        # if both are present, hue takes priority
    colors = []
    color = color if color is not None else sns.color_palette()[0]
    for _hue, batch in zip(hue_order, hue_batches):
        _colors = []
        for _x in x_vals:
            if highlight is not None:
                if _x in highlight:
                    if highlight_color is not None:
                        _colors.append(highlight_color)
                    else:
                        _colors.append(palette.get(_hue, palette.get(_x, color)))
                else:
                    _colors.append(alt_color)
            else:
                _colors.append(palette.get(_hue, palette.get(_x, color)))
        colors.append(_colors)

    # plot kwargs
    if orientation == "horizontal":
        default_plot_kwargs = {"height": 0.8, "linewidth": 1.5, "edgecolor": "w"}
    else:
        default_plot_kwargs = {"width": 0.8, "linewidth": 1.5, "edgecolor": "w"}
    if plot_kwargs is not None:
        default_plot_kwargs.update(plot_kwargs)
    plot_kwargs = default_plot_kwargs

    # legend kwargs
    default_legend_kwargs = {"frameon": True, "loc": "best", "fontsize": 12}
    if legend_kwargs is not None:
        default_legend_kwargs.update(legend_kwargs)
    legend_kwargs = default_legend_kwargs

    # make the plot
    plt.figure(figsize=figsize)
    bottom = np.zeros(len(x_vals))
    for h, d, c in zip(hue_order, batch_data, colors):
        y_vals = np.asarray([d.get(_x, 0) for _x in x_vals])
        if orientation == "horizontal":
            plt.barh(x_vals, y_vals, left=bottom, color=c, label=h, **plot_kwargs)
        else:
            plt.bar(x_vals, y_vals, bottom=bottom, color=c, label=h, **plot_kwargs)
        bottom += y_vals

    # style the plot
    ax = plt.gca()
    if orientation == "horizontal":
        if xlabel is None:
            xlabel = "Frequency (%)" if normalize else "Count"
    else:
        if ylabel is None:
            ylabel = "Frequency (%)" if normalize else "Count"
    ax.set_xlabel(xlabel, fontsize=xlabel_fontsize)
    ax.set_ylabel(ylabel, fontsize=ylabel_fontsize)
    ax.tick_params(
        axis="x", labelsize=xtick_labelsize, labelrotation=xtick_labelrotation
    )
    ax.tick_params(
        axis="y", labelsize=ytick_labelsize, labelrotation=ytick_labelrotation
    )
    for s in ["left", "right", "top"]:
        ax.spines[s].set_visible(False)

    if orientation == "horizontal":
        ax.set_ylim([-0.75, len(x_vals) - 0.25])
    else:
        ax.set_xlim([-0.75, len(x_vals) - 0.25])

    # legend
    if len(hue_batches) > 1 and not hide_legend:
        if orientation == "horizontal":
            leg_kwargs = {"loc": "center left", "bbox_to_anchor": [1.01, 0.5]}
            leg_kwargs.update(legend_kwargs)
        else:
            leg_kwargs = legend_kwargs
        ax.legend(**leg_kwargs)
    if hide_legend or palette is None:
        ax.get_legend().remove()

    # save, show or return the ax
    if figfile is not None:
        plt.tight_layout()
        plt.savefig(figfile)
    elif show:
        plt.show()
    else:
        return ax

