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
import sys
from typing import Union, Iterable, Optional

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl

import seaborn as sns

from natsort import natsorted

from ..utils.color import get_cmap


def donut(
    values,
    counts=None,
    data=None,
    hue=None,
    palette=None,
    color=None,
    cmap=None,
    name=None,
    hue_order=None,
    force_categorical_hue=False,
    lineage_key="lineage",
    figfile=None,
    figsize=(6, 6),
    pairs_only=False,
    alt_color="#F5F5F5",
    edgecolor="white",
    group_singletons=True,
    singleton_color="lightgrey",
    shuffle_colors=False,
    title=None,
    title_fontsize=36,
    title_x=0.5,
    title_y=0.5,
    subtitle=None,
    subtitle_fontsize=20,
    subtitle_x=0.5,
    subtitle_y=0.425,
    show_subtitle=True,
    random_seed=1234,
    width=0.55,
    linewidth=2,
    text_kws={},
    subtext_kws={},
    pie_kws={},
):
    """
    Creates a donut plot of a population of lineages, with arc widths proportional to lineage size.

    .. note::
       For **continuous** hues (for example, AgBC UMI counts), the mean value for each lineage is used. 
       For **boolean** hues (for example, specificity classifications), the lineage is considered ``True`` if
       any lineage member is ``True``. For **categorical** hues (for example, CDR3 length), the most common
       value for each lineage is used. 
    
    Parameters
    ----------
    
    adata : anndata.AnnData  
        Input ``AnnData`` object. ``adata.obs`` must contain a column for the 
        lineage name (`lineage_key`) and, optionally, a `hue` column.
            
    hue : str or dict, optional  
        Can be either the name of a column in ``adata.obs`` or a ``dict`` mapping
        lineage names to hue values. Used to determine the color of each lineage arc. If a ``dict`` 
        is provided, any missing lineage names will still be included in the donut plot but will 
        be colored using `alt_color`. There are three possible classes of hue values:  
             
                - **continuous:** hues that map to a continuous numerical space, identified by all `hue` 
                  values being floating point numbers. An example would be log2-transformed
                  antigen barcode UMI counts. For continuous hues, the mean of all members
                  in a lineage will be plotted.  
                  
                - **boolean:** hues that map to either ``True`` or ``False``. An example would be specificity
                  classification. For boolean hues, if any member of a lineage is ``True``, the
                  entire lineage will be considered ``True``.  
                  
                - **categorical:** hues that map to one of a set of categories. An example would be isotypes. 
                  For categorical hues, the most common value observed in a lineage will 
                  be plotted.  
                  
        Finally, if `hue` is not provided, the lineage name will be considered the `hue`, and
        each lineage will be colored separately.
            
    palette : dict, optional  
        A ``dict`` mapping hue categories to colors. For boolean hue types, if `palette`
        is not provided, `color` will be used for ``True`` and `alt_color` will be used for 
        ``False``. For categorical hue types, if `color` is provided, a monochromatic palette 
        consisting of various shades of `color` will be used. If `color` is not provided,
        ``sns.hls_palette()`` will be used to generate a color palette.  
            
    color : str or list, optional  
        A color name, hex string or RGB list/tuple for coloring the donut plot. For
        boolean hue types, `color` will be used for ``True`` and 'alt_color' will be used for 
        ``False``. For categorical and continuous hue types, a monochromatic palette will 
        be created containing various shades of `color`.  
            
    alt_color : str or list, default='#F5F5F5'    
        A color name, hex string or RGB list/tuple for coloring alternate values 
        (``False`` boolean hues or values not found in `palette`). Default is ``'#F5F5F5'``, which is
        a very light grey.
        
    singleton_color : str or list, default='lightgrey'  
        A color name, hex string or RGB list/tuple for coloring the singleton arc in the donut plot.  
            
    shuffle_colors : bool, default=False  
        If ``True``, colors will be shuffled prior to assignment to hue categories. This
        is primarily useful when the `hue` is the lineage name and a monochromatic palette is used,  
        in order to make it easier to distinguish neighboring arcs on the plot.  
            
    name : str, optional  
        Not currently used.
        
    hue_order : list, optional  
        A list specifying the hue category order. This does not affect the ordering 
        of lineages in the donut plot, just the assignment of colors to `hue` categories. For example,
        when plotting with a monochromatic palette (by providing `color`), `hue_order` will
        order the coloring of `hue` categories from dark to light.
            
    force_categorical_hue : bool, default=False  
        By default, any `hue` categories consisting solely of ``float`` values 
        will be considered continuous and will be colored using a user-supplied colormap (`cmap`) or 
        with a monochromatic color gradient (using `color` as the base color). If ``True``, `hue`
        categories will always be considered categorical.
            
    lineage_key : str, default='lineage'  
        Column in ``adata.obs`` corresponding to the lineage name.  
            
    figfile : str, optional  
        Path to an output figure file. If not provided, the figure will be shown and not saved to file.
            
    figsize : iterable object, default=[6, 6]  
        Figure size, in inches.  
        
    pairs_only : bool, default=False  
        If ``True``, only paired BCR/TCR sequences (containing both heavy/light, alpha/beta or
        delta/gamma chains) will be included.  
            
    edgecolor : str or list, default='white'  
        A color name, hex string or RGB list/tuple for coloring the edges that divide donut arcs.  
            
    random_seed : int, float or str, default=1234  
        Used to set the random seed using ``numpy.random.seed()``. Only applicable when 
        `shuffle_colors` is ``True``, and provided mainly to allow users to recreate plots
        that use shuffled colors (otherwise the shuffle order would be random, thus creating
        a different color order each time the plotting function is called). Default is ``1234``.
            
    width : float, default=0.55  
        Fraction of the donut plot radius that corresponds to the donut 'hole'.  
            
    fontsize : int or float, default=28  
        Fontsize for the sequence count text displayed in the center of the plot.  
            
    linewidth : int or float, default=2  
        Width of the lines separating lineage arcs.  
        
    pie_kws : dict, optional  
        Dictionary containing keyword arguments that will be passed directly to ``ax.pie()``.
            
    text_kws : dict, optional  
        Dictionary containing keyword arguments that will be passed directly to ``ax.text()`` 
        when drawing the text in the center of the plot.
        
    
    """
    # process input data
    _data = {}
    if data is None:
        if counts is None:
            counter = Counter(values)
            values = list(counter.keys())
            counts = list(counter.values())
    else:
        if isinstance(values, str) and values in data.columns:
            values = data[values]
        if counts is None:
            counter = Counter(values)
            values = list(counter.keys())
            counts = list(counter.values())
        elif isinstance(counts, str) and counts in data.columns:
            counts = data[counts]
        if isinstance(hue, str) and hue in data.columns:
            hue = data[hue]
    _data["values"] = values
    _data["counts"] = counts
    if hue is not None:
        _data["hue"] = hue
    df = pd.DataFrame(_data)

    # organize linages into a DataFrame
    ldata = []
    singleton_count = 0
    for i, (l, s) in enumerate(adata.obs[lineage_key].value_counts().items()):
        if s > 1:
            ldata.append({"lineage": l, "size": s, "order": i})
        else:
            singleton_count += 1
    df = pd.DataFrame(ldata)

    # singletons
    singleton_df = pd.DataFrame(
        [
            {
                "lineage": "singletons",
                "size": singleton_count,
                "order": df.shape[0] + 1,
                "hue": "singletons",
                "color": singleton_color,
            },
        ]
    )

    # hue
    if hue is not None:
        if isinstance(hue, dict):
            _hue = [hue.get(l, None) for l in df["lineage"]]
            df["hue"] = _hue
        elif hue in adata.obs:
            if all([isinstance(h, float) for h in adata.obs[hue]]):
                _hue = []
                for l in df["lineage"]:
                    _adata = adata[adata.obs[lineage_key] == l]
                    if _adata:
                        h = np.mean(_adata.obs[hue])
                        _hue.append(h)
                    else:
                        _hue.append(None)
                df["hue"] = _hue
            elif all([isinstance(h, bool) for h in adata.obs[hue]]):
                _hue = []
                for l in df["lineage"]:
                    _adata = adata[adata.obs[lineage_key] == l]
                    if _adata:
                        h = any(_adata.obs[hue])
                        _hue.append(h)
                    else:
                        _hue.append(None)
                df["hue"] = _hue
            else:
                _hue = []
                for l in df["lineage"]:
                    _adata = adata[adata.obs[lineage_key] == l]
                    if _adata:
                        h = _adata.obs[hue].value_counts().index[0]
                        _hue.append(h)
                    else:
                        _hue.append(None)
                df["hue"] = _hue
        else:
            err = "\nERROR: hue must either be the name of a column in adata.obs or a dictionary "
            err += f"mapping lineage names to hue values. You provided {hue}.\n"
            print(err)
            sys.exit()
    else:
        df["hue"] = df["lineage"]
    # set hue type
    if all([isinstance(h, float) for h in df["hue"]]) and not force_categorical_hue:
        hue_type = "continuous"
    elif all([isinstance(h, bool) for h in df["hue"]]):
        hue_type = "boolean"
    elif all([h in df["lineage"] for h in df["hue"]]):
        hue_type = "lineage"
    else:
        hue_type = "categorical"

    # color
    if hue_type == "continuous":
        if cmap is None:
            color = color if color is not None else sns.color_palette()[0]
        cmap = get_cmap(from_color=color) if cmap is None else get_cmap(cmap)
        norm_hue = (df["hue"] - df["hue"].min()) / (df["hue"].max() - df["hue"].min())
        df["color"] = [cmap(nh) for nh in norm_hue]
    elif hue_type == "boolean":
        if palette is None:
            color = color if color is not None else sns.color_palette()[0]
            pos = color
            neg = alt_color
        else:
            pos = palette[True]
            neg = palette[False]
        df["color"] = [pos if h else neg for h in df["hue"]]
    else:
        if palette is None:
            if hue_type == "lineage":
                hue_order = df["hue"]
            else:
                hue_order = (
                    hue_order
                    if hue_order is not None
                    else natsorted(df["hue"].dropna().unique())
                )
            if color is not None:
                colors = _get_monochrome_colors(color, len(hue_order))
                if shuffle_colors:
                    primary = colors[0]
                    secondary = colors[1:]
                    np.random.seed(random_seed)
                    np.random.shuffle(secondary)
                    colors = [primary] + secondary
            else:
                colors = sns.hls_palette(n_colors=len(hue_order))
            palette = {h: c for h, c in zip(hue_order, colors)}
        df["color"] = [palette.get(h, alt_color) for h in df["hue"]]

    # concat the singletons and sort
    df = pd.concat([df, singleton_df], ignore_index=True)
    df = df.sort_values(by="order", ascending=False)

    # make the plot
    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.axis("equal")
    pctdistance = 1 - width / 2
    pie_kwargs = dict(startangle=90, radius=1, pctdistance=1 - width / 2)
    for k, v in pie_kws.items():
        pie_kwargs[k] = v
    pie_kwargs["colors"] = df["color"]
    slices, _ = ax.pie(df["size"], **pie_kwargs)
    plt.setp(slices, width=width, edgecolor=edgecolor)
    for w in slices:
        w.set_linewidth(linewidth)

    # add text to the center of the donut (total sequence count)
    if title is None:
        show_subtitle = False
        title = str(adata.shape[0])
    txt_kwargs = dict(
        size=title_fontsize, color="k", va="center", ha="center", fontweight="bold"
    )
    for k, v in text_kws.items():
        txt_kwargs[k] = v
    ax.text(title_x, title_y, title, transform=ax.transAxes, **txt_kwargs)

    if show_subtitle:
        if subtitle is None:
            subtitle = f"n={adata.shape[0]}"
        subtxt_kwargs = dict(
            size=subtitle_fontsize, color="k", va="center", ha="center"
        )
        for k, v in subtext_kws.items():
            subtxt_kwargs[k] = v
        ax.text(
            subtitle_x, subtitle_y, subtitle, transform=ax.transAxes, **subtxt_kwargs
        )

    plt.tight_layout()

    if figfile is not None:
        plt.savefig(figfile)
    else:
        plt.show()


def _get_monochrome_colors(monochrome_color, n_col):
    cmap = get_cmap(from_color=monochrome_color)
    # this is a bit convoluted, but what's happening is we're getting different colormap
    # values -- which range from 1 (darkest) to 0 (lightest). Calling cmap(i) returns an
    # rgba tuple, but we just need the rbg, so we drop the a. To make sure that one of
    # the colors isn't pure white, we ask np.linspace() for one more value than we need
    # and drop the lightest value
    RGB_tuples = [cmap(i)[:-1] for i in np.linspace(1, 0, n_col + 1)][:-1]
    return RGB_tuples
