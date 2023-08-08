#!/usr/bin/env python
# filename: heatmap.py


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


import sys
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

from .models.data import PlotData
from ..utils.color import get_cmap


def heatmap(
    data: Union[pd.DataFrame, Iterable],
    color: Union[str, Iterable, mpl.colors.Colormap] = "Greys",
    row_colors: Optional[Iterable] = None,
    column_colors: Optional[Iterable] = None,
    scale_color_by_row: bool = False,
    scale_color_by_column: bool = False,
    transform: Optional[str] = None,
    norm: bool = False,
    percent: bool = False,
    norm_axis: Union[int, str] = "columns",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    row_order: Optional[Iterable] = None,
    row_labels: Optional[Iterable] = None,
    row_label_position: str = "left",
    row_label_rotation: float = 0.0,
    row_label_fontsize: float = 14.0,
    show_row_ticks: bool = False,
    column_order: Optional[Iterable] = None,
    column_labels: Optional[Iterable] = None,
    column_label_position: str = "top",
    column_label_rotation: float = 90.0,
    column_label_fontsize: float = 14.0,
    show_column_ticks: bool = False,
    show_values: bool = False,
    values_fmt: Union[str, mpl.ticker.Formatter] = "{x:0.1f}",
    values_fontsize: float = 9.0,
    values_fontweight: str = "normal",
    values_color: Union[str, Iterable, None] = None,
    values_lightcolor: Union[str, Iterable] = "#FFFFFF",
    values_darkcolor: Union[str, Iterable] = "#303030",
    values_horizontalalignment: str = "center",
    values_verticalalignment: str = "center",
    values_color_threshold: Optional[float] = None,
    square: bool = True,
    linewidth: float = 2.0,
    linecolor: str = "#FFFFFF",
    linestyle: str = "-",
    ax: Optional[mpl.axes.Axes] = None,
    figsize: Optional[Iterable] = None,
    figfile: Optional[str] = None,
) -> Optional[mpl.axes.Axes]:
    """
    Creates a heatmap from a pandas DataFrame or an iterable of iterables.

    Parameters
    ----------
    data : Union[pd.DataFrame, Iterable]
        A pandas DataFrame or an iterable of iterables. If a DataFrame is provided,
        the index and columns will be used as the row and column labels, respectively.

    color : Union[str, Iterable, mpl.colors.Colormap], optional
        The color palette to use for the heatmap. Can be any of the following:
          * a ``list`` or ``tuple`` containing RGB or RGBA values
          * a color string, either from `Matplotlib's set of named colors`_ or a hex color value
          * a Matplotlib ``Colormap`` object
          * the name of a registered `Matplotlib Colormap`_
        Default is ``'Greys'``.

    row_colors : Optional[Iterable], optional
        A list of colors to use for each row. The number of colors must match the number of rows.
        Each list item can be anything accepted by `color`. If `row_colors` is provided,
        `color` will be ignored. Default is ``None``.

    column_colors : Optional[Iterable], optional
        A list of colors to use for each column. The number of colors must match the number of
        columns. Each list item can be anything accepted by `color`. If `column_colors` is
        provided, `color` will be ignored. Default is ``None``.

    scale_color_by_row : bool, optional
        If ``True``, the color values will be separately scaled for each row.
        If both `scale_color_by_row` and `scale_color_by_column` are ``False``, the color values
        will be scaled across the entire heatmap. Default is ``False``.

    scale_color_by_column : bool, optional
        If ``True``, the color values will be separately scaled for each column.
        If both `scale_color_by_row` and `scale_color_by_column` are ``False``, the color values
        will be scaled across the entire heatmap. Default is ``False``.

    transform : Optional[str], optional
        The type of transformation to apply to the data. Can be any of the following:
            * ``'log'``: apply a log10 transformation
            * ``'log2'``: apply a log2 transformation
            * ``'log10'``: apply a log10 transformation
            * ``'linear'``: do not apply a transformation
        Default is ``None``, which has the same effect as ``'linear'``.

    vmin : Optional[float], optional
        The minimum value to use for the color scale. If ``None``, the minimum value in the data
        will be used. Default is ``None``.

    vmax : Optional[float], optional
        The maximum value to use for the color scale. If ``None``, the maximum value in the data
        will be used. Default is ``None``.

    row_order : Optional[Iterable], optional
        A list of row labels in the order they should appear in the heatmap. If ``None``, the
        rows will be ordered by the index of the DataFrame. Default is ``None``.

    row_labels : Optional[Iterable], optional
        A list of labels to use for each row. If ``None``, the row labels will be the index of the
        DataFrame. Default is ``None``.

    row_label_position : str, optional
        The position of the row labels. Can be either ``'left'`` or ``'right'``. Default is
        ``'left'``.

    row_label_rotation : float, optional
        The rotation of the row labels in degrees. Default is ``0.0``.

    row_label_fontsize : float, optional
        The font size of the row labels. Default is ``14``.

    show_row_ticks : bool, optional
        If ``True``, ticks will be shown with the row labels. Default is ``False``.

    column_order : Optional[Iterable], optional
        A list of column labels in the order they should appear in the heatmap. If ``None``, the
        columns will be ordered by the columns of the DataFrame. Default is ``None``.

    column_labels : Optional[Iterable], optional
        A list of labels to use for each column. If ``None``, the column labels will be the columns
        of the DataFrame. Default is ``None``.

    column_label_position : str, optional
        The position of the column labels. Can be either ``'top'`` or ``'bottom'``. Default is
        ``'top'``.

    column_label_rotation : float, optional
        The rotation of the column labels in degrees. Default is ``0.0``.

    column_label_fontsize : float, optional
        The font size of the column labels. Default is ``14``.

    show_column_ticks : bool, optional
        If ``True``, ticks will be shown with the column labels. Default is ``False``.

    show_values : bool, optional
        If ``True``, the values will be shown in each cell. Default is ``False``.

    values_fmt : str, optional
        The format string to use for the values. Default is ``'{x:0.2f}'``.

    values_fontsize : float, optional
        The font size of the values. Default is ``8``.

    values_fontweight : str, optional
        The font weight of the values. Default is ``'normal'``.

    values_color : Union[str, Iterable], optional
        The color of the values. Can be a color string, either from
        `Matplotlib's set of named colors`_ or a hex color value, or a list of RGB(A)
        values. Default is ``None``, which will use `values_lightcolor`
        for dark backgrounds and `values_darkcolor` for light backgrounds.

    values_lightcolor : Union[str, Iterable]
        The color of the values when the background is dark. Can be a color string, either from
        `Matplotlib's set of named colors`_ or a hex color value, or a list of RGB(A) values.
        Default is ``'#FFFFFF'`` (white).

    values_darkcolor : Union[str, Iterable]
        The color of the values when the background is light. Can be a color string, either from
        `Matplotlib's set of named colors`_ or a hex color value, or a list of RGB(A) values.
        Default is ``'#303030'`` (dark gray).

    values_horizontalalignment : str, optional
        The horizontal alignment of the values. Can be either ``'left'``, ``'center'``, or
        ``'right'``. Default is ``'center'``.

    values_verticalalignment : str, optional
        The vertical alignment of the values. Can be either ``'top'``, ``'center'``, or
        ``'bottom'``. Default is ``'center'``.

    values_color_threshold : Optional[float], optional
        The threshold to use for determining whether to use `values_lightcolor` or
        `values_darkcolor`. If ``None``, the threshold will be set to the midpoint of the
        background color scale. Default is ``None``.

    square : bool, optional
        If ``True``, heatmap boxes will be square. Default is ``True``.

    linewidth : float, optional
        The width of the lines that will divide each box. Default is ``2.0``.

    linecolor : Union[str, Iterable], optional
        The color of the lines that will divide each box. Can be a color string, either from
        `Matplotlib's set of named colors`_ or a hex color value, or a list of RGB(A) values.
        Default is ``'#FFFFFF'`` (white).

    linestyle : str, optional
        The style of the lines that will divide each box. Default is ``'-'``, which is a
        solid line.

    ax : Optional[mpl.axes.Axes], optional
        The `Matplotlib Axes`_ on which to draw the heatmap. If ``None``, a new figure and axes
        will be created. Default is ``None``.

    figsize : Optional[Tuple[float, float]], optional
        The width and height of the figure in inches. Default is ``None``, which will use
        the dimensions of the input data to infer an appropriate size.

    figfile : Optional[str], optional
        The path to save the figure to. If ``None``, the figure will not be saved and
        the figure ``Axes`` will be returned. Default is ``None``.

    Returns
    -------
    Optional[mpl.axes.Axes]
        The `Matplotlib Axes`_ on which the heatmap was drawn. If ``figfile`` is not ``None``,
        the figure will be saved and the ``Axes`` will be returned.


    .. _Matplotlib's set of named colors:
        https://matplotlib.org/stable/gallery/color/named_colors.html
    -- _Matplotlib Colormap:
        https://matplotlib.org/stable/tutorials/colors/colormaps.html
    .. _Matplotlib Axes:
        https://matplotlib.org/stable/api/axes_api.html

    """
    data = PlotData(data=data, column_labels=column_labels, row_labels=row_labels,)

    # reorder data
    data.reorder(row_order=row_order, column_order=column_order)
    column_labels = data.df.columns.values.tolist()
    row_labels = data.df.index.values.tolist()

    # norm
    if norm or percent:
        data.norm(axis=norm_axis, as_percent=percent)

    # freeze values
    val_df = data.df.copy()
    val_data = val_df.to_numpy()
    n_rows, n_columns = val_df.shape

    # clip data
    data.clip(lower=vmin, upper=vmax)

    # transform
    data.transform(func=transform)

    # scale
    if scale_color_by_row:
        data.scale(axis="rows")
    elif scale_color_by_column:
        data.scale(axis="columns")
    else:
        data.scale()

    # colors
    if row_colors is not None:
        if (n_row_colors := len(row_colors)) < n_rows:
            err = f"\nERROR: the number of row_colors must be equal or greater to the number of rows in the input dataset.\n"
            err += f"The input dataset has {n_rows} rows, but {n_row_colors} colors were provided.\n"
            print(err)
            sys.exit()
        row_cmaps = [get_cmap(c) for c in row_colors]
        color_df = get_color_values(data.df, cmap=row_cmaps, by_row=True)
    elif column_colors is not None:
        if (n_column_colors := len(column_colors)) < n_columns:
            err = f"\nERROR: the number of row_colors must be equal or greater to the number of rows in the input dataset.\n"
            err += f"The input dataset has {n_columns} rows, but {n_column_colors} colors were provided.\n"
            print(err)
            sys.exit()
        column_cmaps = [get_cmap(c) for c in column_colors]
        color_df = get_color_values(data.df, cmap=column_cmaps, by_column=True)
    else:
        cmap = get_cmap(color)
        color_df = get_color_values(data.df, cmap=cmap)

    # reformat the color data
    #
    # color_df.to_numpy() returns a 2D array with tuples (of type "object")
    # as values, rather than a 3D Numpy array. Calling .tolist() on the
    # 2D array retypes the values and allows me to create a 3D
    # array with np.array()
    color_data = np.array(color_df.to_numpy().tolist())

    # set figsize and aspect
    if figsize is None:
        figsize = np.array(color_df.shape[::-1]) / 2
    plt.figure(figsize=figsize)

    # plot
    im = plt.imshow(color_data, aspect="equal" if square else "auto")
    ax = plt.gca()

    # row label positions
    if row_label_position.lower() == "right":
        labelleft = False
        labelright = True
        tickleft = False
        tickright = True if show_row_ticks else False
    else:
        labelleft = True
        labelright = False
        tickleft = True if show_row_ticks else False
        tickright = False
    # row labels
    if row_labels is None:
        row_labels = val_df.index.values
    ax.set_yticks(np.arange(val_df.shape[0]), labels=row_labels)
    ax.tick_params(
        axis="y", labelrotation=row_label_rotation, labelsize=row_label_fontsize,
    )

    # column label positions
    if column_label_position.lower() == "top":
        labeltop = True
        labelbottom = False
        ticktop = True if show_column_ticks else False
        tickbottom = False
    else:
        labeltop = False
        labelbottom = True
        ticktop = False
        tickbottom = True if show_column_ticks else False
    # column labels
    if column_labels is None:
        column_labels = val_df.columns.values
    ax.set_xticks(np.arange(val_df.shape[1]), labels=column_labels)
    ax.tick_params(
        axis="x", labelrotation=column_label_rotation, labelsize=column_label_fontsize,
    )

    # show/hide ticks and labels
    ax.tick_params(
        top=ticktop,
        bottom=tickbottom,
        right=tickright,
        left=tickleft,
        labeltop=labeltop,
        labelbottom=labelbottom,
        labelleft=labelleft,
        labelright=labelright,
    )

    # hide the spines
    ax.spines[:].set_visible(False)

    # draw gridlines separating the heatmap boxes
    ax.set_xticks(np.arange(val_df.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(val_df.shape[0] + 1) - 0.5, minor=True)
    ax.grid(
        which="minor", color=linecolor, linestyle=linestyle, linewidth=linewidth,
    )
    ax.tick_params(which="minor", bottom=False, left=False, top=False, right=False)

    # show values
    if show_values:
        # set the threshold for light/dark text
        if values_color_threshold is not None:
            threshold = values_color_threshold
        else:
            threshold = val_data.max() / 2

        # get the string formatter for the values
        if isinstance(values_fmt, str):
            values_fmt = mpl.ticker.StrMethodFormatter(values_fmt)

        # write the values
        lightcolor = values_color if values_color is not None else values_lightcolor
        darkcolor = values_color if values_color is not None else values_darkcolor
        for i in range(val_data.shape[0]):
            for j in range(val_data.shape[1]):
                c = lightcolor if val_data[i, j] > threshold else darkcolor
                ax.text(
                    x=j,
                    y=i,
                    s=values_fmt(val_data[i, j], None),
                    color=c,
                    horizontalalignment=values_horizontalalignment,
                    verticalalignment=values_verticalalignment,
                    fontsize=values_fontsize,
                    fontweight=values_fontweight,
                )

    # save or return
    if figfile is not None:
        plt.tight_layout()
        plt.savefig(figfile)
    else:
        return ax


def get_color_values(
    df: pd.DataFrame,
    cmap: Union[mpl.colors.Colormap, Iterable[mpl.colors.Colormap]],
    by_row=False,
    by_column=False,
) -> pd.DataFrame:
    """
    Gets colors for each value in a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing values to be colored.

    cmap : Union[mpl.colors.Colormap, Iterable[mpl.colors.Colormap]]
        ``Colormap`` or iterable of ``Colormap`` objects.

    by_row : bool, default=False
        Whether to get colors by row.

    by_column : bool, default=False
        Whether to get colors by column.

    Returns
    -------
    color_df : pd.DataFrame
        DataFrame containing colors for each value in ``df``.
    """
    df = df.copy()
    if by_column:
        return get_colors_by_column(df, cmaps=cmap)
    elif by_row:
        return get_colors_by_row(df, cmaps=cmap)
    else:
        for c in df.columns:
            df[c] = [cmap(v) for v in df[c]]
        return df


def get_colors_by_column(
    df: pd.DataFrame, cmaps: Iterable[mpl.colors.Colormap]
) -> pd.DataFrame:
    """
    Gets colors for each value in a DataFrame, by column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing values to be colored.

    cmaps : Iterable[mpl.colors.Colormap]
        Iterable of ``Colormap`` objects, one for each column in ``df``.

    Returns
    -------
    color_df : pd.DataFrame
        DataFrame containing colors for each value in ``df``.
    """
    for c, cmap in zip(df.columns, cmaps):
        df[c] = [cmap(v) for v in df[c]]
    return df


def get_colors_by_row(
    df: pd.DataFrame, cmaps: Iterable[mpl.colors.Colormap]
) -> pd.DataFrame:
    """
    Gets colors for each value in a DataFrame, by row.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing values to be colored.

    cmaps : Iterable[mpl.colors.Colormap]
        Iterable of ``Colormap`` objects, one for each row in ``df``.

    Returns
    -------
    color_df : pd.DataFrame
        DataFrame containing colors for each value in ``df``.
    """
    df = get_colors_by_column(df.T, cmaps=cmaps)
    return df.T
