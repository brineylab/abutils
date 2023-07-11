#!/usr/bin/env python
# filename: scatter.py


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


import itertools
from typing import Iterable, Union, Optional

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl

import seaborn as sns

from natsort import natsorted

from .data import process_input_data
from .utils import get_inset_axes_bounds
from ..core.sequence import Sequence
from ..utils.color import get_cmap, true_false_palette


__all__ = ["scatter"]


def scatter(
    x: Union[str, Iterable],
    y: Union[str, Iterable],
    hue: Union[str, Iterable, None] = None,
    marker: str = "o",
    data: Optional[pd.DataFrame] = None,
    sequences: Optional[Iterable[Sequence]] = None,
    hue_order: Optional[Iterable] = None,
    force_categorical_hue: bool = False,
    force_continuous_hue: bool = False,
    palette: Union[dict, Iterable, None] = None,
    color: Union[str, Iterable, None] = None,
    cmap: Union[str, mpl.colors.Colormap, None] = None,
    hue_min: Optional[float] = None,
    hue_max: Optional[float] = None,
    under_color: Union[str, Iterable, None] = "whitesmoke",
    size: Union[int, float] = 20,
    alpha: float = 1.0,
    highlight_index: Optional[Iterable] = None,
    highlight_x: Optional[Iterable] = None,
    highlight_y: Optional[Iterable] = None,
    highlight_marker: str = "x",
    highlight_size: Union[int, float] = 90,
    highlight_color: Union[str, Iterable] = "k",
    highlight_name: Optional[str] = None,
    highlight_alpha: float = 0.9,
    plot_kwargs: Optional[dict] = None,
    legend_marker_alpha: Optional[float] = None,
    legend_on_data: bool = False,
    legend_fontsize: Union[int, float] = 12,
    legend_fontweight: str = "bold",
    legend_fontoutline: Union[int, float] = 3,
    legend_label_position_offsets: Optional[dict] = None,
    legend_kwargs: Optional[dict] = None,
    hide_legend: bool = False,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    title_fontsize: Union[int, float] = 20,
    title_fontweight: str = "normal",
    title_loc: str = "center",
    title_pad: Union[int, float, None] = None,
    show_title: bool = False,
    xlabel_fontsize: Union[int, float] = 16,
    ylabel_fontsize: Union[int, float] = 16,
    xtick_labelsize: Union[int, float] = 14,
    ytick_labelsize: Union[int, float] = 14,
    xtick_labelrotation: Union[int, float] = 0,
    ytick_labelrotation: Union[int, float] = 0,
    hide_ticks: bool = False,
    tiny_axis: bool = False,
    tiny_axis_xoffset: Union[int, float, None] = None,
    tiny_axis_yoffset: Union[int, float, None] = None,
    cbar_width: float = 0.35,
    cbar_height: float = 0.05,
    cbar_loc: str = "lower right",
    cbar_orientation: str = "horizontal",
    cbar_bbox_to_anchor: Optional[Iterable] = None,
    cbar_flip_ticks: bool = False,
    cbar_title: Optional[str] = None,
    cbar_title_fontsize: Union[int, float] = 12,
    cbar_title_loc: Optional[str] = None,
    cbar_title_labelpad: float = 8.0,
    hide_cbar: bool = False,
    equal_axes: bool = True,
    ax: Optional[mpl.axes.Axes] = None,
    show: bool = False,
    figsize: Optional[Iterable] = None,
    figfile: Optional[str] = None,
) -> Optional[mpl.axes.Axes]:
    """
    Produces a scatter plot.

    Parameters
    ----------
    x : str or iterable object
        Name of a column in `data` or an iterable of values to be plotted on the
        x-axis. Required.

    y : str or iterable object
        Name of a column in `data` or an iterable of values to be plotted on the
        y-axis. Required.

    hue : str or iterable object, optional
        Name of a column in `data` or an iterable of hue categories.

    marker : str, dict or iterable object, optional
        Marker style for the scatter plot. Accepts any of the following:
          * a `matplotlib marker`_ string
          * a ``dict`` mapping `hue` categories to a `matplotlib marker`_ string
          * a ``list`` of `matplotlib marker`_ strings, which should be the same
              length as `x` and `y`.

    data : pandas.DataFrame, optional
        A ``DataFrame`` object containing the input data. If provided, `x`, `y` and/or `hue` should
        be column names in `data`.

    sequences : iterable of abutils.core.sequence.Sequence, optional
        An iterable of ``Sequence`` objects. If provided, `x`, `y` and/or `hue` should be annotations
        in the ``Sequence`` objects. Alternatively, `x`, `y` and/or `hue` can be an iterable of
        values to be plotted, but must be the same length as `sequences`.

    hue_order : iterable object, optional
        List of `hue` categories in the order they should be plotted. If `hue_order` contains a
        subset of all categories found in `hue`, only the supplied categories will be plotted.
        If not provided, `hue` categories will be plotted in ``natsort.natsorted()`` order.

    force_categorical_hue : bool, default=False
        If ``True``, `hue` data will be treated as categorical, even if the data appear to
        be continuous. This results in `color` being used to color the points rather than `cmap`.

    force_continuous_hue : bool, default=False
        If ``True``, `hue` data will be treated as continuous, even if the data appear to
        be categorical. This results in `cmap` being used to color the points rather than `color`.
        This may produce unexpected results and/or errors if used on non-numerical data.

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
            individual RGB(A) values to a hue category. Wrapping the RGB(A) iterable in a list
            will ensure that the color is interpreted correctly.

        Only used if `hue` contains categorical data (`cmap` is used for continuous data). If not
        provided, the `default Seaborn color palette`_ will be used.

    cmap : str or matplotlib.color.Colormap, default='flare'
        Colormap to be used for continuous `hue` data.

    hue_min : float, default=0
        Minimum value for `hue` when `hue` is continuous. Values below `hue_min` will be set to
        `under_color`.

    hue_max : float, default=1
        Maximum value for `hue` when `hue` is continuous. Values at or above `hue_max` will
        all be colored as the maxmum value in `cmap`.

    under_color : str or list of RGB(A) values
        Separate color for ``0`` values when `hue` is continuous. By default, `cmap` is
        used for all values. An example use would be GEx plots for which visualization is
        improved if ``0`` values are more obviously distinguished from low count values.

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

    legend_on_data : bool, default=False
        Plot legend labels on the data rather than in a separate legend. The X/Y midpoint
        for each legend category is used as the label location.

    legend_fontsize : int or float, default=12
        Fontsize for legend labels.

    legend_fontweight : str, default="normal"
        Fontweight for legend labels. Options are any accepted `Matplotlib text weight`_.

    legend_fontoutline : int or float, default=None
        Width of the outline of legend labels. Only used when `legend_on_data` is ``True``.
        Default is ``None``, which results in no outline.

    legend_marker_alpha : float, default=None
        Opacity for legend markers (or legend labels if `legend_on_data` is ``True``).
        By default, legend markers will use `alpha` and legend labels will be completely
        opaque, equivalent to `legend_marker_alpha` of ``1.0``.

    legend_label_position_offsets : dict, default=None
        A ``dict`` mapping legend labels to ``(x,y)`` coordinates used to offset legend labels.
        Only used when `legend_on_data` is ``True``. Offsets are in relative plot units: ``(0.1, 0.1)``
        would move the label up and to the right by 10% of the overall plot area.

    legend_kwargs : dict, optional
        Dictionary containing keyword arguments that will be passed to ``ax.legend()``.

    hide_legend : bool, default=False
        By default, a plot legend will be shown if multiple hues are plotted. If ``True``,
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

    tiny_axis : bool, default=False
        Plots tiny axis lines in the lower left corner of the plot. Typcally used in
        UMAP plots. If ``True``, ticks and tick labels will be hidden.

    tiny_axis_xoffset : float, default=None
        X-axis offset for the tiny axis.

    tiny_axis_yoffset : float, default=None
        Y-axis offset for the tiny axis.

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

    show : bool, default=False
        If ``True``, plot is shown and the plot ``Axes`` object is not returned. Default
        is ``False``, which does not call ``pyplot.show()`` and returns the ``Axes`` object.

    figsize : iterable object, default=[6, 4]
        List containing the figure size (as ``[x-dimension, y-dimension]``) in inches.

    figfile : str, optional
        Path at which to save the figure file. If not provided, the figure is not saved
        and is either shown (if `show` is ``True``) or the ``Axes`` object is returned.


    Returns
    -------
    ax : mpl.axes.Axes
        If `figfile` is ``None`` and `show` is ``False``, the ``ax`` is returned.
        Otherwise, the return value is ``None``.


    .. _matplotlib marker:
        https://matplotlib.org/stable/api/markers_api.html

    .. _Matplotlib's set of named colors:
        https://matplotlib.org/stable/gallery/color/named_colors.html

    .. _default Seaborn color palette:
        https://seaborn.pydata.org/generated/seaborn.color_palette.html

    .. _Matplotlib text weight
        https://matplotlib.org/stable/tutorials/text/text_props.html

    .. _any valid inset_axes() location:
        https://matplotlib.org/stable/api/_as_gen/mpl_toolkits.axes_grid1.inset_locator.inset_axes.html

    """
    # process input data
    # if data is None:
    #     _data = {}
    #     _data["x"] = x
    #     x = "x"
    #     _data["y"] = y
    #     y = "y"
    #     if hue is not None:
    #         _data["hue"] = hue
    #         hue = "hue"
    #     df = pd.DataFrame(_data)
    # else:
    #     df = data.copy()
    #     if not isinstance(x, str) and len(x) == df.shape[0]:
    #         df["x"] = x
    #         x = "x"
    #     if not isinstance(y, str) and len(y) == df.shape[0]:
    #         df["y"] = y
    #         y = "y"
    #     if hue is not None:
    #         if not isinstance(hue, str) and len(hue) == df.shape[0]:
    #             df["hue"] = hue
    #             hue = "hue"
    df, x, y, hue = process_input_data(
        x=x, y=y, hue=hue, data=data, sequences=sequences
    )

    # figure size
    if figsize is None:
        figsize = [6, 6]

    # hue and color
    continuous_hue = False
    if hue is not None:
        if force_continuous_hue:
            continuous_hue = True
        elif all([isinstance(h, float) for h in df[hue]]) and not force_categorical_hue:
            continuous_hue = True
        if continuous_hue:
            continuous_hue = True
            # set hue min and max values
            if hue_min is None:
                hue_min = np.floor(df[hue].min())
            if hue_max is None:
                hue_max = np.ceil(df[hue].max())
            hue_order = []
            cmap = get_cmap("flare" if cmap is None else cmap)
            normhue = lambda h: (h - hue_min) / (hue_max - hue_min)
            df["color"] = [
                cmap(normhue(h)) if h >= hue_min else under_color for h in df[hue]
            ]
        else:
            if hue_order is None:
                hue_order = natsorted(list(set(df[hue])))
            if palette is not None:
                missing_color = color if color is not None else "lightgrey"
                hue_dict = {h: palette.get(h, missing_color) for h in hue_order}
            else:
                if all([isinstance(h, bool) for h in hue_order]):
                    hue_dict = true_false_palette
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
    # TODO: allow different markers that correspond to marker "categories"
    # much like we can do with hues
    # this might be a bit complex, because pyplot does not currently
    # support marker style assignment by list, so we'd need to make
    # a separate plt.scatter() call for each marker style

    # plot kwargs
    default_plot_kwargs = {"linewidths": 0}
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
            if legend_marker_alpha is not None:
                leg = ax.get_legend()
                for lh in leg.legendHandles:
                    lh.set_alpha(legend_marker_alpha)

    # colorbar
    if continuous_hue and not hide_cbar:
        if cbar_orientation == "horizontal":
            width = max([cbar_width, cbar_height])
            height = min([cbar_width, cbar_height])
        else:
            width = min([cbar_width, cbar_height])
            height = max([cbar_width, cbar_height])
        cbar_bounds = get_inset_axes_bounds(
            cbar_loc, cbar_bbox_to_anchor, width, height
        )
        cbax = ax.inset_axes(cbar_bounds)

        # max_hue = np.ceil(df[hue].max())
        # min_hue = max(0, df[hue].min())
        norm = mpl.colors.Normalize(vmin=hue_min, vmax=hue_max)
        # ticks = [t for t in np.linspace(min_hue, max_hue, num=4)]

        cbar = plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbax,
            orientation=cbar_orientation,
        )
        if cbar_orientation == "horizontal":
            if cbar_title_loc is None:
                cbar_title_loc = "top"
            ticks_position = "bottom" if cbar_flip_ticks else "top"
            cbax.xaxis.set_ticks_position(ticks_position)
            cbax.xaxis.set_label_position(cbar_title_loc)
            cbar.ax.set_xlabel(
                cbar_title,
                fontsize=cbar_title_fontsize,
                labelpad=cbar_title_labelpad,
            )
        else:
            if cbar_title_loc is None:
                cbar_title_loc = "right"
            ticks_position = "left" if cbar_flip_ticks else "right"
            cbax.yaxis.set_ticks_position(ticks_position)
            cbax.yaxis.set_label_position(cbar_title_loc)
            cbar.ax.set_ylabel(
                cbar_title,
                fontsize=cbar_title_fontsize,
                labelpad=cbar_title_labelpad,
            )

    # style the plot
    ax.set_xlabel(xlabel if xlabel is not None else x, fontsize=xlabel_fontsize)
    ax.set_ylabel(ylabel if ylabel is not None else y, fontsize=ylabel_fontsize)
    ax.tick_params(
        axis="x", labelsize=xtick_labelsize, labelrotation=xtick_labelrotation
    )
    ax.tick_params(
        axis="y", labelsize=ytick_labelsize, labelrotation=ytick_labelrotation
    )

    if tiny_axis:
        # get coords for the UMAP-specific axes
        if tiny_axis_xoffset is None:
            tiny_axis_xoffset = -0.05
        if tiny_axis_yoffset is None:
            tiny_axis_yoffset = -0.05
        xmin = df[x].min()
        xmax = df[x].max()
        ymin = df[y].min()
        ymax = df[y].max()
        x_range = abs(xmax - xmin)
        y_range = abs(ymax - ymin)
        x_offset = x_range * tiny_axis_xoffset
        y_offset = y_range * tiny_axis_yoffset
        x_start = xmin + x_offset
        y_start = ymin + y_offset
        x_end = xmin + (x_range / 5) + x_offset
        x_center = x_start + ((x_end - x_start) / 2)
        y_end = ymin + (y_range / 5) + y_offset
        y_center = y_start + ((y_end - y_start) / 2)
        # draw the new "mini" axis lines
        ax.hlines(y_start, x_start, x_end, "k", lw=2)
        ax.vlines(x_start, y_start, y_end, "k", lw=2)
        ax.annotate(
            xlabel,
            xy=(x_center, y_start),
            xytext=(0, -5),
            textcoords="offset points",
            fontsize=xlabel_fontsize,
            ha="center",
            va="top",
        )
        ax.annotate(
            ylabel,
            xy=(x_start, y_center),
            xytext=(-5, 0),
            textcoords="offset points",
            fontsize=ylabel_fontsize,
            rotation="vertical",
            ha="right",
            va="center",
        )
        # remove the normal axis lines
        ax.set_xlabel("", fontsize=0)
        ax.set_ylabel("", fontsize=0)
        for s in ["left", "right", "top", "bottom"]:
            ax.spines[s].set_visible(False)
    else:
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_position(("outward", 10))

    if all([hue is not None, legend_on_data, not continuous_hue]):
        # set up label offsets
        if legend_label_position_offsets is None:
            legend_label_position_offsets = {}
        _xlim = ax.get_xlim()
        _ylim = ax.get_ylim()
        _xrange = _xlim[1] - _xlim[0]
        _yrange = _ylim[1] - _ylim[0]
        # configure the legend font outline
        if legend_fontoutline is not None:
            path_effect = [
                mpl.patheffects.withStroke(linewidth=legend_fontoutline, foreground="w")
            ]
        else:
            path_effect = None
        # add the on-data legend
        for h in df[hue].unique():
            _df = df[df[hue] == h]
            xoffset, yoffset = legend_label_position_offsets.get(h, (0, 0))
            hue_x = _df[x].median() + (xoffset * _xrange)
            hue_y = _df[y].median() + (yoffset * _yrange)
            ax.text(
                hue_x,
                hue_y,
                h,
                c="k",
                alpha=legend_marker_alpha,
                weight=legend_fontweight,
                verticalalignment="center",
                horizontalalignment="center",
                fontsize=legend_fontsize,
                path_effects=path_effect,
            )
        # hide the "real" legend
        hide_legend = True

    if hide_legend:
        l = ax.get_legend()
        if l is not None:
            l.remove()

    if hide_ticks or tiny_axis:
        ax.set_xticks([])
        ax.set_yticks([])

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


# def get_inset_axes_bounds(loc, bbox_to_anchor, width, height):
#     if bbox_to_anchor is None:
#         loc_dict = {
#             "upper left": [0, 1 - height],
#             "upper center": [0.5 - width / 2, 1 - height],
#             "upper right": [1 - width, 1 - height],
#             "center left": [0, 0.5 - height / 2],
#             "center": [0.5 - width / 2, 0.5 - height / 2],
#             "center right": [1 - width, 0.5 - height / 2],
#             "lower left": [0, 0],
#             "lower center": [0.5 - width / 2, 0],
#             "lower right": [1 - width, 0],
#         }
#         x0, y0 = loc_dict.get(loc, [0, 0])
#     else:
#         x, y = bbox_to_anchor[:2]
#         loc_dict = {
#             "upper left": [x, y - height],
#             "upper center": [x - width / 2, y - height],
#             "upper right": [x - width, y - height],
#             "center left": [x, y - height / 2],
#             "center": [x - width / 2, y - height / 2],
#             "center right": [x - width, y - height / 2],
#             "lower left": [x, y],
#             "lower center": [x - width / 2, y],
#             "lower right": [x - width, y],
#         }
#         x0, y0 = loc_dict.get(loc, [0, 0])
#     return [x0, y0, width, height]
