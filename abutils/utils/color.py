#!/usr/bin/env python
# filename: color.py


#
# Copyright (c) 2015 Bryan Briney
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


from typing import Union, Tuple, Optional, Iterable

import seaborn as sns

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import cm, colors
from matplotlib.colors import Colormap, LinearSegmentedColormap, to_rgba_array


cmaps = {"heatmap": sns.diverging_palette(240, 10, as_cmap=True)}


def get_cmap(
    c: Union[Colormap, str, Iterable],
    dark: bool = False,
    zero_color: Union[str, Iterable, None] = None,
    n: int = 256,
    minval: float = 0.0,
    maxval: float = 1.0,
    name: Optional[str] = None,
) -> Colormap:
    """
    Gets a matplotlib ``Colormap``.

    Parameters
    ----------
    c : ``str``, ``Colormap`` or ``tuple``
        Can be one of several things:
            * matplotlib ``Colormap``
            * hex code
            * `Matplotlib color`_
            * RGB tuple
        If a single color is provided, a ``Colormap`` will be built from 
        white to the provided color (or from the color to black if `dark` 
        is ``True``). 
    
    dark : bool, default=False
        If ``True``, build a ``Colormap`` from the provided color from the color 
        to black. By default, the ``Colormap`` will be built from the provided
        color to white.

    zero_color : any Matplotlib color, default=None
        Alternate color for the zero value of the ``Colormap``. Produces results
        similar to the default colormap in 10x Genomics' Loupe browser, in which 
        the zero count value for GEX plots is light grey and less easily confused with 
        near-zero counts.

    minval : float, default=0.0
        Truncate the lower bound of the ``Colormap``. Must be a ``float`` less than 1.0,
        which represents the fraction of the ``Colormap`` at which to truncate.

    maxval : float, default=1.0
        Truncate the upper bound of the ``Colormap``. Must be a ``float`` greater than 
        zero and less than or equal to 1.0, which represents the fraction of the 
        ``Colormap`` to truncate.

    n : int, default=256
        Number of steps in the ``Colormap``.

    name : str, default=None
        Name of the ``Colormap``.

    
    Returns
    -------
    cmap : ``Colormap``


    .. _Matplotlib color
        https://matplotlib.org/stable/tutorials/colors/colors.html
    """
    if isinstance(c, Colormap):
        cmap = c
    elif isinstance(c, str) and c in cm.cmap_d:
        cmap = plt.get_cmap(c)
    elif dark:
        cmap = sns.dark_palette(c, as_cmap=True)
    else:
        cmap = sns.light_palette(c, as_cmap=True)
    if zero_color is not None:
        cmap = cmap(np.linspace(0.1, 1, n))
        colors = [to_rgba_array(zero_color)] + list(cmap)
        cmap = LinearSegmentedColormap.from_list(
            name if name is not None else "", colors
        )
    cmap = LinearSegmentedColormap.from_list(name, cmap(np.linspace(minval, maxval, n)))
    return cmap


def monochrome_palette(
    color: Union[str, Iterable], n_colors: int = 10, include_white: bool = False
) -> list:
    """
    Returns a monochromatic palette of colors, from `color` to white.
    
    Parameters
    ----------
    color : str or Iterable
        Color from which the monochromatic palette will be created. 
        Can be one of several things:
            * hex code
            * `Matplotlib color`_
            * RGB tuple 

    n_colors : int, default=10
        Number of colors in the palette.

    include_white : bool, default=False
        Whether or not to include white in the palette. White is not 
        included by default, meaning all colors in the palette will be 
        a shade of `color`.


    Returns
    -------
    palette : list
        A ``list`` of RGBA ``tuple``s


    .. _Matplotlib color
        https://matplotlib.org/stable/tutorials/colors/colors.html
    """
    cmap = get_cmap(color)
    if include_white:
        colors = [cmap(i) for i in np.linspace(1, 0, n_colors)]
    else:
        colors = [cmap(i) for i in np.linspace(1, 0, n_colors + 1)][:-1]
    return colors


def cmap_from_color(color: Union[str, Tuple], dark: bool = False) -> Colormap:
    """
    Generates a matplotlib colormap from a single color. Colormap will be built, 
    by default, from white to ``color``.

    Parameters
    ----------
    color : str or tuple
        Can be one of several things:
            1. Hex code
            2. HTML color name
            3. RGB tuple

    dark : bool, default=False
        If ``True``, colormap will be built from ``color`` to
        black. Default is ``False``, which builds a colormap from
        white to ``color``.

    Returns
    -------
        colormap : ``Colormap``
    """
    if dark:
        return sns.dark_palette(color, as_cmap=True)
    else:
        return sns.light_palette(color, as_cmap=True)


def hex_to_rgb(hex_string):
    rgb = colors.hex2color(hex_string)
    return tuple([int(255 * x) for x in rgb])


def rgb_to_hex(rgb_tuple):
    div = 1 if all([v <= 1.0 for v in rgb_tuple]) else 255
    return colors.rgb2hex([1.0 * x / div for x in rgb_tuple])


def hls(n_colors, hue=0.01, lightness=0.6, saturation=0.65):
    return sns.hls_palette(n_colors, h=hue, l=lightness, s=saturation)


def husl(n_colors, hue=0.01, saturation=0.9, lightness=0.65):
    return sns.husl_palette(n_colors, h=hue, s=saturation, l=lightness)


def truncate_colormap(
    cmap: Colormap, minval: float = 0.0, maxval: float = 1.0, n: int = 256
) -> Colormap:
    """
    Truncates a colormap, such that the new colormap consists of
    ``cmap[minval:maxval]``.

    If maxval is larger than minval, the truncated colormap will be reversed.

    Args:

    	cmap : colormap): Colormap to be truncated.

    	minval (float): Lower bound. Should be a float betwee 0 and 1.

    	maxval (float): Upper bound. Should be a float between 0 and 1

    	n (int): Number of colormap steps. Default is ``256``.

    Returns:

    	colormap: A matplotlib colormap

    http://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
    """
    cmap = get_cmap(cmap)
    name = "%s-trunc-%.2g-%.2g" % (cmap.name, minval, maxval)
    return colors.LinearSegmentedColormap.from_list(
        name, cmap(np.linspace(minval, maxval, n))
    )


# def stack_colormap(lower, upper, n=256):
#     """
#     Stacks two colormaps (``lower`` and ``upper``) such that
#     low half -> ``lower`` colors, high half -> ``upper`` colors

#     Args:

#     	lower (colormap): colormap for the lower half of the stacked colormap.

#     	upper (colormap): colormap for the upper half of the stacked colormap.

#     	n (int): Number of colormap steps. Default is ``256``.
#     """
#     A = get_cmap(lower)
#     B = get_cmap(upper)
#     name = "%s-%s" % (A.name, B.name)
#     lin = np.linspace(0, 1, n)
#     return array_cmap(np.vstack((A(lin), B(lin))), name, n=n)


# def get_cmap(cmap=None, name=None, from_color=None, dark=False, n=256):
#     # """
#     # Generates a matplotlib colormap.

#     # cmap can be one of several things:
#     #     - a name ('Blues', 'BuGn_r') of a built-in colormap
#     #     - a cmap
#     #     - a filename, np.loadtxt() n x 3 or 4  ints 0..255 or floats 0..1
#     #     - a numpy array
#     # See http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps or in IPython, plt.cm.<tab>

#     # An optional name for the colormap can be provided (name), as well as the number
#     # of cmap steps (n).

#     # Alternatively, to make a colormap using a single color, provide the color
#     # via from_color. By default, the supplied color will be the dark end of
#     # the cmap, with white as the lightest color. To reverse (use the input
#     # color as the lighest value and black as the darkest), set dark=True.

#     # Returns a matplotlib colormap object.
#     # """
#     if from_color is not None:
#         return cmap_from_color(from_color, dark)
#     elif cmap is None:
#         err = "You must provide either cmap or from_color"
#         raise RuntimeError(err)
#     if isinstance(cmap, colors.Colormap):
#         return cmap
#     if isinstance(cmap, str):
#         if cmap in cm.cmap_d:
#             return plt.get_cmap(cmap)  # "Blues" ...
#         A = np.loadtxt(cmap, delimiter=None)  # None: white space
#         name = name or cmap.split("/")[-1].split(".")[0]  # .../xx.csv -> xx
#     else:
#         A = cmap  # numpy array or array-like
#     return array_cmap(A, name, n=n)


# def array_cmap(A, name=None, n=256):
#     """ numpy array -> a cmap, matplotlib.colors.Colormap
#         n x 3 or 4  ints 0 .. 255 or floats 0 ..1
#     """
#     A = np.asanyarray(A)
#     assert A.ndim == 2 and A.shape[1] in (
#         3,
#         4,
#     ), "array must be n x 3 or 4, not %s" % str(A.shape)
#     Amin, Amax = A.min(), A.max()
#     if A.dtype.kind == "i":
#         assert 0 <= Amin < Amax <= 255, "Amin %d  Amax %d must be in 0 .. 255" % (
#             Amin,
#             Amax,
#         )
#         A = A / 255.0  # not /=
#     else:
#         assert 0 <= Amin < Amax <= 1, "Amin %g  Amax %g must be in 0 .. 1" % (
#             Amin,
#             Amax,
#         )
#     return colors.LinearSegmentedColormap.from_list(name or "noname", A, N=n)

