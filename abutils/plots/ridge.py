#!/usr/bin/env python
# filename: ridge.py


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


import warnings

warnings.filterwarnings("ignore")

import itertools
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from natsort import natsorted


def ridge(
    categories: Union[str, Iterable, None] = None,
    values: Union[str, Iterable, None] = None,
    data: Optional[pd.DataFrame] = None,
    order: Optional[Iterable] = None,
    palette: Union[dict, Iterable, None] = None,
    alt_color: Union[Iterable, str] = "lightgrey",
    alpha: float = 1.0,
    linewidth: float = 0.0,
    outlinewidth: float = 1.5,
    xaxis_linewidth: float = 1.0,
    xaxis_linecolor: Union[Iterable, str, None] = None,
    xlabel: str = "UMI count ($\mathregular{log_2}$)",
    xlabel_fontsize: int = 12,
    ylabel_fontsize: int = 11,
    bw: Union[str, float] = "scott",
    adjust_hspace: float = 0.1,
    category_label_fontweight: str = "bold",
    category_label_xoffset: float = 0.25,
    xmin: Union[int, float, None] = None,
    xmax: Union[int, float, None] = None,
    aspect: float = 15.0,
    height: float = 0.5,
    figfile: Optional[str] = None,
    show: bool = True,
):
    """
    Draws a ridge plot.

    Parameters
    ----------
    categories : Union[str, Iterable, None], optional
        Category information. Can be any of the following:
          - if `data` is provided and `values` is not provided, `categories`
            must be a list containing one or more column names in `data`
          - if both `data` and `values` are also provided, `categories` should be a
            column name in `data`
          - if `data` is not provided, `categories` must be a list of category names

    values : Union[str, Iterable, None], optional
        Value information. Can be any of the following:
          - if `data` and `categories` are provided, `values` should be a column name in `data`
          - if `data` is not provided, `values` must be a list of values
        Note that `values` cannot be provided without `categories`.

    data : Optional[pd.DataFrame], optional
        Input ``DataFrame``. If not provided, ``categories`` and ``values`` must be provided.

    order : Optional[Iterable], optional
        Order of categories in the plot. If not provided, categories will be sorted using ``natsort``.

    palette : Union[dict, Iterable, None], optional
        Color palette. Can be any of the following:
            - a dictionary mapping categories to colors
            - a list of colors
        If a ``dict`` is provided, any categories not included in the dictionary will be assigned
        `alt_color`. If a ``list`` is provided, the colors will be assigned to the categories in order,
        with colors reused if there are more categories than colors. If not provided, a default palette
        (created using ``sns.hls_palette()``) will be used.

    alt_color : Union[Iterable, str], optional
        Color to use for categories not included in the palette. If a ``list`` is provided, the color

    alpha : float, optional
        Alpha value for the density curves. Default is 1.0.

    linewidth : float, optional
        Line width for the data line on the density curves. Default is 0.0.

    outlinewidth : float, optional
        Line width for the white outline of the density curves. The purpose is to provide visual
        separation if/when the density curves of adjacent ridge plots overlap. Default is 1.5.

    xaxis_linewidth : float, optional
        Line width for the x-axis line. Default is 1.0.

    xaxis_linecolor : Union[Iterable, str, None], optional
        Color for the x-axis line. If not provided, the color will be the same as the color
        used to plot the density curve.

    xlabel : str, optional
        Label for the x-axis. Default is "UMI count ($\mathregular{log_2}$)".

    xlabel_fontsize : int, optional
        Font size for the x-axis label. Default is 12.

    ylabel_fontsize : int, optional
        Font size for the y-axis labels. Default is 11.

    bw : Union[str, float], optional
        Bandwidth for the density curves. Can be any of the following:
            - a float
            - the name of a reference rule (see ``scipy.stats.gaussian_kde``)
        Default is "scott".

    adjust_hspace : float, optional
        Adjust the vertical space between subplots. Default is 0.1.

    category_label_fontweight : str, optional
        Font weight for the category labels. Default is "bold".

    category_label_xoffset : float, optional
        Horizontal offset for the category labels, as a fraction of the total plot width.
        Default is 0.25.

    xmin : Union[int, float, None], optional
        Minimum value for the x-axis. If not provided, the minimum value will be the data minimum.

    xmax : Union[int, float, None], optional
        Maximum value for the x-axis. If not provided, the maximum value will be the data maximum.

    aspect : float, optional
        Aspect ratio for the plot. Default is 15.0.

    height : float, optional
        Height of each individual ridge plot, in inches. Default is 0.5.

    figfile : Optional[str], optional
        If provided, the figure will be saved to this file. Default is None.

    show : bool, optional
        If True, the figure will be displayed. Default is True.

    Returns
    -------
    sns.FacetGrid
        If `show` is ``False`` and `figfile` is not provided, the ``FacetGrid`` object
        will be returned.

    """
    # process input data
    if data is None:
        # if just `categories` and `values` are provided, we create a dataframe
        d = {}
        d["category"] = categories
        d["values"] = values
        df = pd.DataFrame(d)
        categories = "category"
        values = "values"
    elif all([data is not None, categories is not None, values is None]):
        # if `data` is provided without `values`, `categories` must be a list of column names
        # we melt the data to get the dataframe into the correct format
        df = data.copy()
        if not all([c in df.columns for c in categories]):
            raise ValueError(
                f"Invalid column names for categories: {categories}. "
                f"If `data` is provided without `values`, `categories` must be a list of column names. "
                f"Available columns: {df.columns}"
            )
        df = df[categories]
        df = df.melt(value_vars=categories, var_name="category", value_name="values")
        categories = "category"
        values = "values"
    else:
        # if `data` is provided with `values` and `categories`, we check to see if they are column names
        # if not (and if the length matches the number of rows in `data`), we create new columns
        df = data.copy()
        if not isinstance(categories, str) and len(categories) == df.shape[0]:
            df["category"] = categories
            categories = "category"
        if not isinstance(values, str) and len(values) == df.shape[0]:
            df["values"] = values
            values = "values"

    # categories and order
    if order is None:
        order = natsorted(df[categories].unique())
    n_categories = len(order)
    df = df.query(f"{categories} in @order")

    # colors
    if palette is None:
        palette = {
            cat: color
            for cat, color in zip(order, sns.hls_palette(n_categories, l=0.5, s=0.8))
        }
    elif palette is not None and isinstance(palette, (list, tuple)):
        iter_colors = itertools.cycle(palette)
        palette = {cat: color for cat, color in zip(order, iter_colors)}
    elif palette is not None and isinstance(palette, dict):
        palette = {cat: palette.get(cat, alt_color) for cat in order}
    else:
        raise ValueError(
            f"Invalid palette: {palette}. Must be an iterable of colors or dict of category: color pairs."
        )
    df["color"] = df[categories].map(palette)

    # make the plot
    g = sns.FacetGrid(
        df,
        row=categories,
        hue=categories,
        aspect=aspect,
        height=height,
        palette=palette,
        row_order=order,
        hue_order=order,
    )
    g.map(
        sns.kdeplot, values, clip_on=False, shade=True, alpha=alpha, lw=linewidth, bw=bw
    )
    g.map(sns.kdeplot, values, clip_on=False, color="w", lw=outlinewidth, bw=bw)
    g.map(plt.axhline, y=0, lw=xaxis_linewidth, color=xaxis_linecolor, clip_on=False)

    # add labels to the individual ridge plots
    def label(x, color, label):
        ax = plt.gca()
        ax.text(
            0,
            0.3,
            label,
            fontweight=category_label_fontweight,
            color=color,
            fontsize=ylabel_fontsize,
            ha="left",
            va="center",
            transform=ax.transAxes,
        )

    g.map(label, categories)

    # style the plot
    if xmin is None:
        xmin = df[values].min()
    if xmax is None:
        xmax = df[values].max()
    category_label_xoffset = (xmax - xmin) * category_label_xoffset
    g.fig.subplots_adjust(hspace=adjust_hspace)
    g.set_titles("")
    g.set(xticks=range(int(np.floor(xmin,)), int(np.ceil(xmax)) + 1, 2,))
    g.set(xlim=(xmin - category_label_xoffset, xmax + 0.25))
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    # set and locate the X axis label
    g.set(ylabel="")
    g.set(xlabel=xlabel)
    xlabel_position = (((xmax - xmin) / 2) + category_label_xoffset) / (
        (xmax - xmin) + category_label_xoffset
    )
    for ax in g.axes.flat:
        ax.set_xlabel(ax.get_xlabel(), x=xlabel_position, fontsize=xlabel_fontsize)

    # save or show the plot
    if figfile:
        g.savefig(figfile)
    elif show:
        plt.show()
    return g
