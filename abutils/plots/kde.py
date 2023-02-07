#!/usr/bin/env python
# filename: kde.py


#
# Copyright (c) 2022 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to de al in the Software without restriction,
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

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import seaborn as sns

from natsort import natsorted

from abutils.utils.color import get_cmap




def kde(
    x=None,
    y=None, 
    hue=None, 
    marker="o",
    data=None,
    hue_order=None, 
    only_scatter_hue=False,
    force_categorical_hue=False,
    palette=None,
    color=None,
    cmap=None,
    size=20,
    alpha=0.6,
    highlight_index=None,
    highlight_x=None,
    highlight_y=None,
    highlight_marker="x",
    highlight_size=90,
    highlight_color="k",
    highlight_name=None,
    highlight_alpha=0.9,  
    kde_kwargs=None,
    scatter_kwargs=None,
    legend_kwargs=None,
    hide_legend=False,
    xlabel=None,
    ylabel=None,
    title=None,
    title_fontsize=20,
    title_fontweight='normal',
    title_loc='center',
    title_pad=None,
    show_title=False,
    xlabel_fontsize=16,
    ylabel_fontsize=16,
    xtick_labelsize=14,
    ytick_labelsize=14,
    xtick_labelrotation=0,
    ytick_labelrotation=0,
    cbar_width=0.35,
    cbar_height=0.05,
    cbar_loc="lower right",
    cbar_orientation="horizontal",
    cbar_bbox_to_anchor=None,
    cbar_flip_ticks=False,
    cbar_title=None,
    cbar_title_fontsize=12,
    hide_cbar=False,
    equal_axes=True,
    ax=None,
    show=False,
    figsize=None,
    figfile=None,
)




    # process input data
    if data is None:
        _data = {}
        _data['x'] = x
        x = 'x'
        _data['y'] = y
        y = 'y'
        if hue is not None:
            _data['hue'] = hue
            hue = 'hue'
        df = pd.DataFrame(_data)
    else:
        df = data.copy()
    
    # figure size
    if figsize is None:
        figsize = [6, 6]
        
    # hue and color
    continuous_hue = False
    if hue is not None:
        if all([isinstance(h, float) for h in df[hue]]) and not force_categorical_hue:
            continuous_hue = True
            hue_order = []
            if cmap is None:
                cmap = sns.color_palette("flare", as_cmap=True)
            else:
                cmap = plt.get_cmap(cmap)
            min_hue = max(0, df[hue].min())
            max_hue = np.ceil(df[hue].max())
            df["color"] = [cmap((h - min_hue) / (max_hue - min_hue)) for h in df[hue]]
        else:
            if hue_order is None:
                hue_order = natsorted(list(set(df[hue])))
            if palette is not None:
                missing_color = color if color is not None else 'lightgrey'
                hue_dict = {h: palette.get(h, missing_color) for h in hue_order}
            else:
                n_colors = max(1, len(hue_order))
                if color is None:
                    color = sns.color_palette(n_colors=n_colors)
                if len(color) < n_colors:
                    color = itertools.cycle(color)
                hue_dict = {h: c for h, c in zip(hue_order, color)}
            df["color"] = [hue_dict[h] for h in df[hue]]
    else:
        hue_order = []
        if isinstance(color, str) and color in df.columns:
            df["color"] = df[color]
        elif color is not None:
            df["color"] = [color] * df.shape[0]
        else:
            df["color"] = [sns.color_palette()[0]] * df.shape[0]










def scatter(
    x=None,
    y=None, 
    hue=None, 
    marker="o",
    data=None,
    hue_order=None, 
    force_categorical_hue=False,
    palette=None,
    color=None,
    cmap=None,
    size=20,
    alpha=0.6,
    highlight_index=None,
    highlight_x=None,
    highlight_y=None,
    highlight_marker="x",
    highlight_size=90,
    highlight_color="k",
    highlight_name=None,
    highlight_alpha=0.9,  
    plot_kwargs=None,
    legend_kwargs=None,
    hide_legend=False,
    xlabel=None,
    ylabel=None,
    title=None,
    title_fontsize=20,
    title_fontweight='normal',
    title_loc='center',
    title_pad=None,
    show_title=False,
    xlabel_fontsize=16,
    ylabel_fontsize=16,
    xtick_labelsize=14,
    ytick_labelsize=14,
    xtick_labelrotation=0,
    ytick_labelrotation=0,
    cbar_width=0.35,
    cbar_height=0.05,
    cbar_loc="lower right",
    cbar_orientation="horizontal",
    cbar_bbox_to_anchor=None,
    cbar_flip_ticks=False,
    cbar_title=None,
    cbar_title_fontsize=12,
    hide_cbar=False,
    equal_axes=True,
    ax=None,
    show=False,
    figsize=None,
    figfile=None,
):
    '''
    Produces a scatter plot.

    Parameters
    ----------

    x : str or iterable object  
        Name of a column in `data` or an iterable of values to be plotted on the 
        x-axis. Required.

    y : str or iterable object  
        Name of a column in `data` or an iterable of values to be plotted on the 
        y-axis..

    hue : str or iterable object, optional  
        Name of a column in `data` or an iterable of hue categories to be used to 
        group data into stacked bars. If not provided, an un-stacked bar plot is created.  
        
    marker : str, dict or iterable object, optional  
        Marker style for the scatter plot. Accepts any of the following:
          * a `matplotlib marker`_ string
          * a ``dict`` mapping `hue` categories to a `matplotlib marker`_ string
          * a ``list`` of `matplotlib marker`_ strings, which should be the same 
              length as `x` and `y`.

    data : pandas.DataFrame, optional  
        A ``DataFrame`` object containing the input data. If provided, `x` and/or `y` should 
        be column names in `data`.   
        
    hue_order : iterable object, optional  
        List of `hue` categories in the order they should be plotted. If `hue_order` contains a 
        subset of all categories found in `hue`, only the supplied categories will be plotted.  
        If not provided, `hue` categories will be plotted in ``natsort.natsorted()`` order.  
        
    force_categorical_hue : bool, default=False  
        If ``True``, `hue` data will be treated as categorical, even if the data appear to 
        be continuous. This results in `color` being used to color the points rather than `cmap`.   
        
    palette : dict, optional  
        Dictionary mapping `hue`, `x` or `y` names to colors. If if keys in `palette` match 
        more than one category, `hue` categories take priority. If `palette` is not provided, 
        bars are colored using `color` (if `hue` is ``None``) or a palette is generated 
        automatically using ``sns.hls_palette()``.  

    color : str or iterable object, optional. 
        Color for the plot markers. Can be any of the following:
          * a ``list`` or ``tuple`` containing RGB or RGBA values
          * a color string, either from `Matplotlib's set of named colors`_ or a hex color value
          * the name of a column in `data` that contains color values
          * a ``list`` of colors (either as strings or RGB tuples) to be used for `hue` categories. 
            If `colors` is shorter than the list of hue categories, colors will be reused.  
        
        .. tip::
            If a single RGB or RGBA ``list`` or ``tuple`` is provided and `hue` is also supplied, 
            there may be unexpected results as ``scatter()`` will attempt to map each of the 
            individual RGB(A) values to a hue category.   
        
        Only used if `hue` contains categorical data (`cmap` is used for continuous data). If not 
        provided, the `default Seaborn color palette`_ will be used. 
        
    cmap : str or matplotlib.color.Colormap, default='flare'   
        Colormap to be used for continuous `hue` data.  
        
    size : str or float or iterable object, default=20  
        Size of the scatter points. Either a   

    alpha : float, default=0.6  
        Alpha of the scatter points.  
        
    highlight_index : iterable object, optional  
        An iterable of index names (present in `data.index`) of points to be highlighted on 
        the plot. If provided, `highlight_x` and `highlight_y` are ignored.

    highlight_x : iterable object, optional  
        An iterable of x-values for highlighted points. Also requires `highlight_y`.
        
    highlight_y : iterable object, optional  
        An iterable of y-values for highlighted points. Also requires `highlight_x`.
        
    highlight_marker : str, default='x'  
        Marker style to be used for highlight points. Accepts any `matplotlib marker`_. 

    highlight_size : int, default=90  
        Size of the highlight marker.

    highlight_color : string or list of color values, default='k'
        Color of the highlight points.

    highlight_name : str, optional  
        Name of the highlights, to be used in the plot legend. If not supplied,
        highlight points will not be included in the legend.
        
    highlight_alpha : float, default=0.9  
        Alpha of the highlight points. 

    plot_kwargs : dict, optional  
        Dictionary containing keyword arguments that will be passed to ``pyplot.scatter()``.

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

    cbar_width : int, default=35  
        Width of the colorbar. Only used for continuous `hue` types.  

    cbar_height : int, default=5  
        Height of the colorbar. Only used for continuous `hue` types.  

    cbar_loc : str or iterable object, default='lower right'  
        Location of the colorbar. Accepts `any valid inset_axes() location`_.

    cbar_orientation : str, default='horizontal'  
        Orientation of the colorbar. Options are ``'horizontal'`` and ``'vertical'``.

    cbar_bbox_to_anchor : list or tuple, optional
        bbox_to_anchor for the colorbar. Used in combination with `cbar_loc` to provide 
        fine-grained positioning of the colorbar.

    cbar_flip_ticks : bool, default=False  
        Flips the position of colorbar ticks. Ticks default to the bottom if `cbar_orientation` 
        is  ``'horizontal'`` and the left if  `cbar_orientation` is ``'vertical'``.  

    cbar_title : str, optional  
        Colorbar title. If not provided, `hue` is used.  

    cbar_title_fontsize : int or float, default=12  
        Fontsize for the colorbar title.  
        
    hide_cbar : bool, default=False. 
        If ``True``, the color bar will be hidden on plots with continuous `hue` values.  
        
    equal_axes : bool, default=True
        If ```True```, the the limits of the x- and y-axis will be equal.

    ax : mpl.axes.Axes, default=None
        Pre-existing axes for the plot. If not provided, a new axes will be created.

    show :bool, default=False  
        If ``True``, plot is shown and the plot ``Axes`` object is not returned. Default
        is ``False``, which does not call ``pyplot.show()`` and returns the ``Axes`` object.

    figsize : iterable object, default=[6, 4]  
        List containing the figure size (as ``[x-dimension, y-dimension]``) in inches. 

    figfile : str, optional  
        Path at which to save the figure file. If not provided, the figure is not saved
        and is either shown (if `show` is ``True``) or the ``Axes`` object is returned.  
        
        
    .. _matplotlib marker: 
        https://matplotlib.org/stable/api/markers_api.html  

    .. _Matplotlib's set of named colors:
        https://matplotlib.org/stable/gallery/color/named_colors.html
        
    .. _default Seaborn color palette: 
        https://seaborn.pydata.org/generated/seaborn.color_palette.html  

    .. _any valid inset_axes() location: 
        https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.axes_grid1.inset_locator.inset_axes.html
        
    '''
    # process input data
    if data is None:
        _data = {}
        _data['x'] = x
        x = 'x'
        _data['y'] = y
        y = 'y'
        if hue is not None:
            _data['hue'] = hue
            hue = 'hue'
        df = pd.DataFrame(_data)
    else:
        df = data.copy()
    
    # figure size
    if figsize is None:
        figsize = [6, 6]
        
    # hue and color
    continuous_hue = False
    if hue is not None:
        if all([isinstance(h, float) for h in df[hue]]) and not force_categorical_hue:
            continuous_hue = True
            hue_order = []
            if cmap is None:
                cmap = sns.color_palette("flare", as_cmap=True)
            else:
                cmap = plt.get_cmap(cmap)
            min_hue = max(0, df[hue].min())
            max_hue = np.ceil(df[hue].max())
            df["color"] = [cmap((h - min_hue) / (max_hue - min_hue)) for h in df[hue]]
        else:
            if hue_order is None:
                hue_order = natsorted(list(set(df[hue])))
            if palette is not None:
                missing_color = color if color is not None else 'lightgrey'
                hue_dict = {h: palette.get(h, missing_color) for h in hue_order}
            else:
                n_colors = max(1, len(hue_order))
                if color is None:
                    color = sns.color_palette(n_colors=n_colors)
                if len(color) < n_colors:
                    color = itertools.cycle(color)
                hue_dict = {h: c for h, c in zip(hue_order, color)}
            df["color"] = [hue_dict[h] for h in df[hue]]
    else:
        hue_order = []
        if isinstance(color, str) and color in df.columns:
            df["color"] = df[color]
        elif color is not None:
            df["color"] = [color] * df.shape[0]
        else:
            df["color"] = [sns.color_palette()[0]] * df.shape[0]
            
    # markers
    # TODO
    
            
    # plot kwargs
    default_plot_kwargs = {'linewidths': 0}
    if plot_kwargs is not None:
        default_plot_kwargs.update(plot_kwargs)
    plot_kwargs = default_plot_kwargs
    
    
    # scatterplot
    if ax is None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
    if hue_order:
        for h in hue_order[::-1]:
            d = df[df[hue] == h]
            ax.scatter(
                d[x],
                d[y],
                c=d["color"],
                s=size,
                marker=marker,
                alpha=alpha,
                label=h,
                **plot_kwargs,
            )
    else:
        ax.scatter(
            df[x],
            df[y],
            c=df["color"],
            s=size,
            marker=marker,
            alpha=alpha,
            **plot_kwargs,
        )

    # highlighted points
    highlight = any(
        [
            highlight_index is not None,
            all([highlight_x is not None, highlight_y is not None]),
        ]
    )
    if highlight:
        if highlight_index is not None:
            hi_index = [h for h in highlight_index if h in df.index.values]
            hi_df = df.loc[hi_index]
            highlight_x = hi_df[x]
            highlight_y = hi_df[y]
        plt.scatter(
            highlight_x,
            highlight_y,
            zorder=10,
            s=highlight_size,
            c=highlight_color,
            alpha=highlight_alpha,
            marker=highlight_marker,
            label=highlight_name,
        )
    
    # legend
    if not continuous_hue:
        if hue is not None and not hide_legend:
            default_legend_kwargs = {"frameon": True, "loc": "best", "fontsize": 12}
            if legend_kwargs is not None:
                default_legend_kwargs.update(legend_kwargs)
            legend_kwargs = default_legend_kwargs
            ax.legend(**legend_kwargs)
    
    # colorbar
    elif not hide_cbar:
        if cbar_orientation == 'horizontal':
            width = max([cbar_width, cbar_height])
            height = min([cbar_width, cbar_height])
        else:
            width = min([cbar_width, cbar_height])
            height = max([cbar_width, cbar_height])
        cbar_bounds = get_inset_axes_bounds(
            cbar_loc,
            cbar_bbox_to_anchor,
            width,
            height
        )
        cbax = ax.inset_axes(cbar_bounds)

        max_hue = np.ceil(df[hue].max())
        min_hue = max(0, df[hue].min())
        norm = mpl.colors.Normalize(vmin=min_hue, vmax=max_hue)
        # ticks = [t for t in np.linspace(min_hue, max_hue, num=4)]
        
        cbar = plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbax,
            orientation=cbar_orientation,
        )
        if cbar_orientation == "horizontal":
            ticks_position = "bottom" if cbar_flip_ticks else "top"
            cbax.xaxis.set_ticks_position(ticks_position)
            cbax.xaxis.set_label_position(ticks_position)
            cbar.ax.set_xlabel(
                cbar_title,
                fontsize=cbar_title_fontsize,
            )
        else:
            ticks_position = "left" if cbar_flip_ticks else "right"
            cbax.yaxis.set_ticks_position(ticks_position)
            cbax.yaxis.set_label_position(ticks_position)
            cbar.ax.set_ylabel(
                cbar_title,
                fontsize=cbar_title_fontsize,
            )

    # style the plot
    ax.set_xlabel(xlabel if xlabel is not None else x, fontsize=xlabel_fontsize)
    ax.set_ylabel(ylabel if ylabel is not None else y, fontsize=ylabel_fontsize)
    ax.tick_params(
        axis="x", labelsize=xtick_labelsize, labelrotation=xtick_labelrotation
    )
    ax.tick_params(axis="y", labelsize=ytick_labelsize, labelrotation=ytick_labelrotation)

    for spine in ["right", "top"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_position(("outward", 10))

    if equal_axes:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        axlim = [min([xlim[0], ylim[0]]), max([xlim[1], ylim[1]])]
        ax.set_xlim(axlim)
        ax.set_ylim(axlim)
        
    if show_title and title is not None:
        ax.set_title(
            title,
            loc=title_loc,
            pad=title_pad,
            fontsize=title_fontsize,
            fontweight=title_fontweight,
        )

    # save, show or return the ax
    if figfile is not None:
        plt.tight_layout()
        plt.savefig(figfile)
    elif show:
        plt.show()
    else:
        return ax
















def feature_kde(
    data,
    x,
    y,
    hue=None,
    hue_order=None,
    colors=None,
    thresh=0.1,
    show_scatter=True,
    scatter_size=5,
    scatter_alpha=0.2,
    fill=False,
    kde_fill_alpha=0.7,
    kde_line_alpha=1.0,
    highlight_index=None,
    highlight_x=None,
    highlight_y=None,
    highlight_marker="x",
    highlight_size=90,
    highlight_color="k",
    highlight_name=None,
    highlight_alpha=0.8,
    xlabel=None,
    ylabel=None,
    equal_axes=True,
    legend_kwargs=None,
    return_ax=False,
    figsize=[6, 6],
    figfile=None,
    **kwargs,
):
    """
    Produces a 2-dimensional KDE plot of two features.

    Parameters
    ----------
    data : anndata.AnnData or pandas.DataFrame  
        An ``AnnData`` object or a ``DataFrame`` containing the input data. Required.

    x : str
        Name of the column in `data` containing the feature to be plotted on the x-axis. Required.

    y : str  
        Name of the column in `data` containing the feature to be plotted on the y-axis. Required.

    hue : str, optional  
        Name of the column in `data` containing categories for hue values. ``hue```` categories 
        will each be plotted as differently colored densities on the same plot. 

    hue_order : iterable object, optional  
        Iterable of hue categories, in the order they should be plotted and listed
        in the legend. If `hue_order` contains only a subset of the categories
        present in ``data[hue]`` or ``data.obs[hue]``, only the categories supplied in `hue_order`
        will be plotted.

    colors : iterable object, optional  
        List of colors to be used for `hue` categories. If `colors` is shorter than the 
        list of hue categories, colors will be reused. If not provided, the 
        `default Seaborn color palette`_ will be used. 

    thresh : float, default=0.1  
        Threshold for the KDE, as a fraction of the overall dataset.
        
    show_scatter : bool, default=True  
        Show a scatterplot beneath the transparent KDE plot.

    scatter_size : int or float, default=5  
        Size of the scatter points.

    scatter_alpha : float, default=0.2  
        Alpha of the scatter points.

    fill : bool, default=True  
        Whether or not to fill the KDE KDE plot. If ``False``, only the KDE boundary lines
        will be plotted.

    kde_fill_alpha : float, default=0.7  
        Alpha for the filled KDE plot. Ignored if `fill` is ``False``.  
        
    kde_line_alpha : float, default=1.0  
        Alpha for the KDE boundary lines.

    highlight_index : iterable object, optional  
        An iterable of index names (present in `data`) of points to be highlighted on 
        the KDE plot. If provided, `highlight_x` and `highlight_y` are ignored.

    highlight_x : iterable object, optional  
        An iterable of x-values for highlighted points. Also requires `highlight_y`.
        
    highlight_y : iterable object, optional  
        An iterable of y-values for highlighted points. Also requires `highlight_x`.
        
    highlight_marker : str, default='x'  
        Marker style to be used for highlight points. Accepts any `matplotlib marker`_. 

    highlight_size : int, default=90  
        Size of the highlight marker.

    highlight_color : string or list of color values, default='k'
        Color of the highlight points.

    highlight_name : str, optional  
        Name of the highlights, to be used in the plot legend. If not supplied,
        highlight points will not be included in the legend.
        
    highlight_alpha : float, default=0.8  
        Alpha of the highlight points.

    xlabel : str, optional  
        Label for the x-axis. By default, the value for `x` is used.

    ylabel : str, optional  
        Label for the y-axis. By default, the value for `y` is used.

    equal_axes : bool, default=True
        If ```True```, the the limits of the x- and y-axis will be equal.
        
    legend_kwargs : dict, optional  
        Dictionary of legend keyword arguments, which will be passed to ``ax.legend()``.

    return_ax : bool, default=False  
        If ``True``, return the plot's ``ax`` object. Will not show or save the plot.

    figsize : list, default=[6, 6]  
        A list containg the dimensions of the plot, in inches.

    figfile : str, optional  
        Path to which the figure will be saved. If not provided, the figure will be 
        shown but not saved to file.

    kwargs  
        All other keyword arguments are passed to ``seaborn.kdeplot()``.

    
    .. _default Seaborn color palette: 
        https://seaborn.pydata.org/generated/seaborn.color_palette.html  

    .. _matplotlib marker: 
        https://matplotlib.org/stable/api/markers_api.html

    """

    # input data
    if isinstance(data, AnnData):
        _data = {}
        for var in [x, y, hue]:
            if var is not None:
                if any([var in data.obs.columns.values, var in data.var_names]):
                    _data[var] = data.obs_vector(var)
                else:
                    print(
                        '"{}" was not found in the supplied AnnData object.'.format(var)
                    )
                    return
        df = pd.DataFrame(_data, index=data.obs_names)
    else:
        _data = {}
        for var in [x, y, hue]:
            if var is not None:
                if var in data.columns.values:
                    _data[var] = data[var]
                else:
                    print('"{}" is not a column in the supplied dataframe'.format(x))
                    return
        df = pd.DataFrame(_data, index=data.index.values)

    # hue
    if hue is not None:
        if hue_order is None:
            hue_order = natsorted(list(set(df[hue])))
        df = df[df[hue].isin(hue_order)]
    else:
        hue_order = []

    # colors
    n_colors = max(1, len(hue_order))
    if colors is None:
        colors = sns.hls_palette(n_colors=n_colors)

    plt.figure(figsize=figsize)

    # scatterplots
    if show_scatter:
        if hue_order:
            for h, c in zip(hue_order, colors):
                d = df[df[hue] == h]
                plt.scatter(
                    d[x], d[y], c=[c], s=scatter_size, alpha=scatter_alpha, linewidths=0
                )
        else:
            plt.scatter(
                df[x],
                df[y],
                c=[colors[0]],
                s=scatter_size,
                alpha=scatter_alpha,
                linewidths=0,
            )

    # kdeplot
    if fill:
        if hue_order:
            sns.kdeplot(
                data=df,
                x=x,
                y=y,
                hue=hue,
                fill=True,
                alpha=kde_fill_alpha,
                hue_order=hue_order,
                palette=colors,
                thresh=thresh,
                **kwargs,
            )
        else:
            sns.kdeplot(
                data=df,
                x=x,
                y=y,
                fill=True,
                alpha=kde_fill_alpha,
                color=colors[0],
                thresh=thresh,
                **kwargs,
            )
    if hue_order:
        ax = sns.kdeplot(
            data=df,
            x=x,
            y=y,
            hue=hue,
            alpha=kde_line_alpha,
            hue_order=hue_order,
            palette=colors,
            thresh=thresh,
            **kwargs,
        )
    else:
        ax = sns.kdeplot(
            data=df,
            x=x,
            y=y,
            alpha=kde_line_alpha,
            color=colors[0],
            thresh=thresh,
            **kwargs,
        )

    # highlighted points
    highlight = any(
        [
            highlight_index is not None,
            all([highlight_x is not None, highlight_y is not None]),
        ]
    )
    if highlight:
        if highlight_index is not None:
            hi_index = [h for h in highlight_index if h in df.index.values]
            hidata = df.loc[hi_index]
            highlight_x = hidata[x]
            highlight_y = hidata[y]
        plt.scatter(
            highlight_x,
            highlight_y,
            zorder=10,
            s=highlight_size,
            c=highlight_color,
            alpha=highlight_alpha,
            marker=highlight_marker,
        )

    # legend
    legend_params = {"loc": "best", "title": None, "fontsize": 12, "frameon": False}
    legend_params.update(legend_kwargs if legend_kwargs is not None else {})
    legend_labels = hue_order
    if fill:
        handles = []
        for c in colors:
            f = Patch(fc=c, alpha=kde_fill_alpha / 3)
            e = Patch(ec=c, fill=False, lw=1.5)
            handles.append((f, e))
    else:
        handles = [Line2D([0], [0], color=c) for c in colors]
    if highlight_name is not None:
        legend_labels.append(highlight_name)
        handles.append(
            Line2D(
                [0],
                [0],
                marker=highlight_marker,
                color="w",
                mec=highlight_color,
                mfc=highlight_color,
                ms=highlight_size / 10,
            )
        )
    ax.legend(handles, legend_labels, **legend_params)

    # style the plot
    ax.set_xlabel(xlabel if xlabel is not None else x, fontsize=16)
    ax.set_ylabel(ylabel if ylabel is not None else y, fontsize=16)
    ax.tick_params(axis="both", labelsize=13)

    for spine in ["right", "top"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_position(("outward", 10))

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