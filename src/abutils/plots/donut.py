#!/usr/bin/env python
# filename: donut.py


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
from typing import Union, Iterable, Optional, Callable

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl

import seaborn as sns

from natsort import natsorted, natsort_keygen

from ..utils.color import get_cmap, monochrome_palette


def donut(
    values: Union[Iterable, str],
    counts: Union[Iterable, str, None] = None,
    hue: Union[Iterable, str, None] = None,
    data: Optional[pd.DataFrame] = None,
    palette: Optional[dict] = None,
    color: Union[str, Iterable, None] = None,
    cmap: Optional[mpl.colors.Colormap] = None,
    order: Optional[Iterable] = None,
    hue_order: Optional[Iterable] = None,
    hue_agg: Union[str, Callable] = "most common",
    sort_by: str = "count",
    sort_key: Callable = natsort_keygen(),
    sort_descending: bool = False,
    force_categorical_hue: bool = False,
    force_continuous_hue: bool = False,
    alt_color: Union[str, Iterable] = "#F5F5F5",
    edgecolor: Union[str, Iterable] = "white",
    singleton_color: Union[str, Iterable, None] = None,
    group_singletons: bool = True,
    shuffle_colors: bool = False,
    random_seed: Union[int, float, str] = 1234,
    title: Optional[str] = None,
    title_fontsize: Union[int, float] = 36,
    title_x: Union[int, float] = 0.5,
    title_y: Union[int, float] = 0.5,
    subtitle: Optional[str] = None,
    subtitle_fontsize: Union[int, float] = 20,
    subtitle_x: Union[int, float] = 0.5,
    subtitle_y: Union[int, float] = 0.425,
    show_subtitle: bool = True,
    width: Union[int, float] = 0.55,
    linewidth: Union[int, float] = 2,
    plot_kwargs: Optional[dict] = None,
    text_kwargs: Optional[dict] = None,
    subtext_kwargs: Optional[dict] = None,
    ax: Optional[mpl.axes.Axes] = None,
    show: bool = False,
    figsize: Optional[Iterable] = None,
    figfile: Optional[str] = None,
) -> Optional[mpl.axes.Axes]:
    """
    Creates a donut plot of a population of lineages, with arc widths proportional to lineage size.

    .. note::
       For **continuous** hues (for example, AgBC UMI counts), the mean value for each lineage is used. 
       For **boolean** hues (for example, specificity classifications), the lineage is considered ``True`` if
       any lineage member is ``True``. For **categorical** hues (for example, CDR3 length), the most common
       value for each lineage is used. 
    
    Parameters
    ----------
    values : str or Iterable
        Name of a column in `data` or an iterable of values to be used as segments of 
        the donut plot. Required.

    counts : str or Iterable
        Name of a column in `data` or an iterable of counts used to determine the size
        of each donut segment. If not provided, counts will be computed by calling 
        ``Counter()`` on `values`.
            
    hue : str or dict, optional  
        Can be either the name of a column in `data` or a ``dict`` mapping
        values to hues. If a ``dict`` is provided, any missing values will 
        still be included in the donut plot but will be colored using `alt_color`. 
        There are three possible classes of hues:  
             
                - **continuous:** hues that map to a continuous numerical space. Hues are
                  automatically calssified as continuous if all `hue` values are ``float``s. 
                  If `force_continuous_hue` is ``True``, the hue will be considered 
                  continuous regardless of whether all hue values are ``float``s.   
                  
                - **boolean:** hues that map to either ``True`` or ``False``. For boolean 
                  hues, if the `hue` for any member of a group of values is ``True``, 
                  the entire group of values will be considered ``True``.  
                  
                - **categorical:** hues that map to one of a set of categories. For 
                  categorical hues, the most common hue for a group of values will be 
                  used for the entire group.   
                  
        Finally, if `hue` is not provided, `values` will be considered the `hue`, and
        each unique value will be colored separately.
            
    palette : dict, default=None  
        A ``dict`` mapping `hue` categories to colors. For boolean `hue` types, if `palette`
        is not provided, `color` will be used for ``True`` and `alt_color` will be used for 
        ``False``. For categorical `hue` types, if neither `color` nor `palette` is provided, 
        ``seaborn.hls_palette`` will be used to generate a color palette. For continuous `hue` 
        types, a ``Colormap`` will be used and `palette` is ignored.  
            
    color : str or list, default=None  
        A `Matplotlib color`_ name, hex string or RGB list/tuple for coloring the donut plot. For
        boolean `hue` types, `color` will be used for ``True`` and 'alt_color' will be used for 
        ``False``. For categorical and continuous hue types, a monochromatic palette will 
        be created that contains various shades from `color` to white (but not including white). 
        If not provided, the first color from ``seaborn.color_palette`` is used.  

    cmap : str or matplotlib.color.Colormap, default=None   
        Colormap to be used for continuous `hue` data. If not provided, `color` is used to 
        make a monochromatic ``Colormap``. If `color` is not provided, Matplotlib's built-in
        ``"flare"`` ``Colormap`` is used.

    order : list, default=None  
        A list specifying the order of donut segments. Elements in `order` should be `values` 
        categories. If not provided, segments will be ordered using `sort_by` and `sort_key`.

    hue_order : list, default=None  
        A list specifying the order of hue categories. This does not affect the ordering 
        of segments in the donut plot, just the assignment of colors to `hue` categories. 
        For example, when plotting with a monochromatic palette (by providing `color`), 
        `hue` categories will be colored from dark to light in `hue_order` order. To re-order 
        donut segments by hue categories, use `sort_by`.

    hue_agg : str or callable, default="most common"
        Method by which to aggregate `hue` values. Can be a ``str`` or a function. Built-in
        options are:
            - **most common:** selects the most common hue value. Suitable for all hue types.
            - **max**: selects the maximum value. Suitable primarily for numerical hue types.
            - **min**: selects the minimum value. Suitable primarily for numerical hue types.
            - **mean**: selects the average value. Suitable primarily for numerical hue types.
            - **median**: selects the median value. Suitable primarily for numerical hue types.
        If a callable is provided, it should accept a list of hue values and return a single 
        value that will be used for the entire category.

    sort_by : str, default="count"
        Field to be used to sort the donut segments. Default is ``"count"``, which sorts 
        the segments by size. Options are:
            - **count:** sorts the segments by the `counts` category (which is the segment size)
            - **value:** sorts the segments by the `values` category 
            _ **hue:** sorts the segments by the `hue` category

    sort_key : callable, default=natsort.natsort_keygen()
        Callable function that determines the sorting method. Will be passed directly to 
        ``sorted()`` as the ``key`` argument.  

    sort_descending : bool, default=False
        If ``True``, sort donut segments in descending order.  

    force_categorical_hue : bool, default=False  
        If ``True``, `hue` categories will always be considered categorical. 

    force_continuous_hue : bool, default=False  
        If ``True``, `hue` categories will always be considered continuous. 
            
    alt_color : str or list, default='#F5F5F5'    
        A color name, hex string or RGB list/tuple for coloring alternate values 
        (``False`` boolean hues or values not found in `palette`). Default is 
        ``'#F5F5F5'``, which maps to a very light grey.  

    edgecolor : str or list, default='white'  
        A `Matplotlib color`_ name, hex string or RGB list/tuple for coloring the edges that 
        separate donut segments. 
        
    singleton_color : str or list, default=None  
        A `Matplotlib color`_ name, hex string or RGB list/tuple for coloring the singleton 
        segment of the donut plot. Only used if `group_singletons` is ``True``. If not provided,
        `alt_color` is used.  

    group_singletons : bool, default=False
        Group all singletons (`values` categories with a `count` of ``1``) into a single, 
        undivided segment.
            
    shuffle_colors : bool, default=False  
        If ``True``, colors will be shuffled prior to assignment to hue categories. This
        can be useful when a monochromatic palette is used, in order to make it easier to 
        distinguish neighboring segments on the plot.   
            
    random_seed : int, float or str, default=1234  
        Used to set the random seed using ``numpy.random.seed()``. Only applicable when 
        `shuffle_colors` is ``True``, and provided mainly to allow users to recreate plots
        that use shuffled colors (otherwise the shuffle order would be random, thus creating
        a different color order each time the plotting function is called).  

    title : str, default=None
        Title, printed in the center of the donut "hole". If not provided, the total 
        count of all elements in `values` is used.

    title_fontsize : float or int, default=36
        Fontsize for the title.

    title_x : float, default=0.5
        X-axis position for the center of the `title` text, in relative coordinates. 
        Should be in the range of 0-1. Default is ``0.5``, which is centered horizontally.

    title_y : float, default=0.5
        Y-axis position for the center of the `title` text, in relative coordinates. 
        Should be in the range of 0-1. Default is ``0.5``, which is centered vertically.  

    subtitle : str, default=None
        Title, printed in the center of the donut "hole". If not provided, the total 
        count of all elements in `values` is used, formatted like ``"n={count}"``.

    subtitle_fontsize : float or int, default=20
        Fontsize for the title.

    subtitle_x : float, default=0.5
        X-axis position for the center of the `subtitle` text, in relative coordinates. 
        Should be in the range of 0-1. Default is ``0.5``, which is centered horizontally.

    subtitle_y : float, default=0.425
        Y-axis position for the center of the `title` text, in relative coordinates. 
        Should be in the range of 0-1. Default is ``0.425`` which, when paired with the 
        default `title_fontsize` and `title` coordinates, will be positioned just under
        the `title`.  

    show_subtitle : bool, default=True
        Whether to show `subtitle`. Only used if `title` is provided, since the default 
        data use to populate `title` and `subtitle` is the same.
             
    width : float, default=0.55  
        Fraction of the donut plot radius that corresponds to the donut 'hole'.  
            
    linewidth : int or float, default=2  
        Width of the lines separating donut segments.  
        
    plot_kwargs : dict, default=None  
        Dictionary containing keyword arguments that will be passed directly to ``ax.pie()``.
            
    text_kwargs : dict, default=None  
        Dictionary containing keyword arguments that will be passed directly to ``ax.text()`` 
        when drawing the title text in the center of the plot.

    subtext_kwargs : dict, default=None  
        Dictionary containing keyword arguments that will be passed directly to ``ax.text()`` 
        when drawing the subtitle text in the center of the plot.

    ax : mpl.axes.Axes, default=None
        Pre-existing axes for the plot. If not provided, a new axes will be created.

    show : bool, default=False  
        If ``True``, plot is shown and the plot ``Axes`` object is not returned. Default
        is ``False``, which does not call ``pyplot.show()`` and returns the ``Axes`` object.

    figsize : iterable object, default=[6, 6]  
        Figure size, in inches. 
            
    figfile : str, optional  
        Path to an output figure file. If not provided, the figure will be shown and not saved to file. 
    

    Returns
    -------
    ax : mpl.axes.Axes
        If `figfile` is ``None`` and `show` is ``False``, the ``ax`` is returned. 
        Otherwise, the return value is ``None``. 
        
    
    .. _Matplotlib color
        https://matplotlib.org/stable/tutorials/colors/colors.html

    """
    # process input data
    d = {}
    # if data is None, then counts/values must be iterables
    if data is None:
        if counts is None:
            counter = Counter(values)
            d["value"] = list(counter.keys())
            d["count"] = list(counter.values())
    # if data isn't None, then we need to determine whether each of
    # counts/values/hue is in data is an iterable, or isn't provided
    else:
        if isinstance(values, str) and values in data.columns:
            values = data[values]
            d["value"] = values
        if counts is None:
            counter = Counter(values)
            d["value"] = list(counter.keys())
            d["count"] = list(counter.values())
        elif isinstance(counts, str) and counts in data.columns:
            d["count"] = data[counts]

    # hue aggregation functions
    most_common_agg_func = lambda x: Counter(x).most_common()[0][0]
    hue_agg_funcs = {
        "max": max,
        "min": min,
        "mean": np.mean,
        "median": np.median,
        "most_common": most_common_agg_func,
    }

    # process hues
    if hue is not None:
        if isinstance(hue, str) and hue in data.columns:
            hue = data[hue]
            hue_vals = {}
            for v, h in zip(values, hue):
                if v not in hue_vals:
                    hue_vals[v] = []
                hue_vals[v].append(h)
            if all([isinstance(h, bool) for h in hue]):
                hue_agg_func = any
            else:
                hue_agg_func = hue_agg_func.get(hue_agg, most_common_agg_func)
            agg_hue_vals = []
            for v in d["value"]:
                h = hue_vals[v]
                if len(h) == 1:
                    agg_hue_vals.append(h[0])
                else:
                    agg_hue_vals.append(hue_agg_func(h))
            d["hue"] = agg_hue_vals
        else:
            d["hue"] = [hue.get(v, None) for v in d["value"]]
    else:
        d["hue"] = d["value"]
    df = pd.DataFrame(d)

    # singletons
    if group_singletons:
        total_count = df.shape[0]
        df = df[df["count"] > 1]
        singleton_count = total_count - df.shape[0]
        singleton_color if singleton_color is not None else alt_color
        singleton_df = pd.DataFrame(
            [
                {
                    "value": "singletons",
                    "count": singleton_count,
                    "order": df.shape[0] + 1,
                    "hue": "singletons",
                    "color": singleton_color,
                },
            ]
        )

    # set the plotting order
    if order is not None:
        order_lookup = {v: i for i, v in enumerate(order)}
        df["order"] = [order_lookup[v] for v in df["value"]]
    else:
        if sort_by in ["hue", "hues"]:
            sort_by = "hue"
        elif sort_by in ["value", "values"]:
            sort_by = "value"
        else:
            sort_by = "count"
        order = sorted(list(set(df[sort_by])), key=sort_key, reverse=sort_descending)
        order_lookup = {v: i for i, v in enumerate(order)}
        df["order"] = [order_lookup[v] for v in df[sort_by]]

    # set hue type
    if force_continuous_hue:
        hue_type = "continuous"
    elif all([isinstance(h, float) for h in df["hue"]]) and not force_categorical_hue:
        hue_type = "continuous"
    elif all([isinstance(h, bool) for h in df["hue"]]):
        hue_type = "boolean"
    else:
        hue_type = "categorical"

    # color
    if hue_type == "continuous":
        if cmap is None:
            if color is not None:
                cmap = get_cmap(color)
            else:
                cmap = get_cmap("flare")
        else:
            cmap = get_cmap(cmap)
        norm_hue = (df["hue"] - df["hue"].min()) / (df["hue"].max() - df["hue"].min())
        df["color"] = [cmap(nh) for nh in norm_hue]
    elif hue_type == "boolean":
        if palette is None:
            pos = color if color is not None else sns.color_palette()[0]
            neg = alt_color
        else:
            pos = palette[True]
            neg = palette[False]
        df["color"] = [pos if h else neg for h in df["hue"]]
    else:
        if palette is None:
            if hue_order is None:
                hue_order = natsorted(df["hue"].dropna().unique())
            if color is not None:
                colors = monochrome_palette(color, n_colors=len(hue_order))
                if shuffle_colors:
                    primary = colors[0]
                    secondary = colors[1:]
                    np.random.seed(random_seed)
                    np.random.shuffle(secondary)
                    colors = [primary] + secondary
            else:
                colors = sns.hls_palette(n_colors=len(hue_order))
                if shuffle_colors:
                    np.random.seed(random_seed)
                    np.random_shuffle(colors)
            palette = {h: c for h, c in zip(hue_order, colors)}
        df["color"] = [palette.get(h, alt_color) for h in df["hue"]]

    # concat the singletons and sort
    if group_singletons:
        df = pd.concat([df, singleton_df], ignore_index=True)
    df = df.sort_values(by="order", ascending=False)

    # make the plot
    if ax is None:
        plt.figure(figsize=figsize if figsize is not None else [6, 6])
        ax = plt.gca()
    ax.axis("equal")
    plot_kws = {"startangle": 90, "radius": 1, "pctdistance": 1 - width / 2}
    if plot_kwargs is not None:
        for k, v in plot_kwargs.items():
            plot_kws[k] = v
    plot_kws["colors"] = df["color"]
    slices, _ = ax.pie(df["count"], **plot_kws)
    plt.setp(slices, width=width, edgecolor=edgecolor)
    for w in slices:
        w.set_linewidth(linewidth)

    # add text to the center of the donut (default is total number of items)
    if title is None:
        show_subtitle = False
        title = str(len(values))
    txt_kws = {
        "size": title_fontsize,
        "color": "k",
        "va": "center",
        "ha": "center",
        "fontweight": "bold",
    }
    if text_kwargs is not None:
        for k, v in text_kwargs.items():
            txt_kws[k] = v
    ax.text(title_x, title_y, title, transform=ax.transAxes, **txt_kws)

    # add the subtitle (default is the total number of items)
    if show_subtitle:
        if subtitle is None:
            subtitle = f"n={len(values)}"
        subtxt_kws = {
            "size": subtitle_fontsize,
            "color": "k",
            "va": "center",
            "ha": "center",
        }
        if subtext_kwargs is not None:
            for k, v in subtext_kwargs.items():
                subtxt_kws[k] = v
        ax.text(subtitle_x, subtitle_y, subtitle, transform=ax.transAxes, **subtxt_kws)

    plt.tight_layout()
    if figfile is not None:
        plt.savefig(figfile)
    elif show:
        plt.show()
    else:
        return ax


# def _get_monochrome_colors(monochrome_color, n_col):
#     cmap = get_cmap(monochrome_color)
#     # this is a bit convoluted, but what's happening is we're getting different colormap
#     # values -- which range from 1 (darkest) to 0 (lightest). Calling cmap(i) returns an
#     # rgba tuple, but we just need the rbg, so we drop the a. To make sure that one of
#     # the colors isn't pure white, we ask np.linspace() for one more value than we need
#     # and drop the lightest value
#     RGB_tuples = [cmap(i)[:-1] for i in np.linspace(1, 0, n_col + 1)][:-1]
#     return RGB_tuples
