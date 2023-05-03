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
import matplotlib as mpl
import matplotlib.pyplot as plt

from natsort import natsorted

import pycircos


def ring(
    values: Union[str, dict, Iterable],
    categories: Optional[Union[str, Iterable]] = None,
    data: Optional[pd.DataFrame] = None,
    color: Optional[Union[str, Iterable]] = None,
    palette: Optional[Union[str, Iterable, dict]] = None,
    raxis_range: Optional[tuple] = None,
    linewidth: Optional[float] = 1.0,
    facecolor: Optional[str] = "#FFFFFF",
    label_visible: Optional[bool] = False,
    arc_start: Optional[float] = 0,
    arc_end: Optional[float] = 360,
    title: Optional[str] = None,
    title_fontsize: Optional[float] = 36,
    title_fontweight: Optional[str] = "heavy",
    subtitle: Optional[str] = None,
    subtitle_fontsize: Optional[float] = 20,
    subtitle_fontweight: Optional[str] = "normal",
    show_category_labels: Optional[bool] = True,
    category_label_fontsize: Optional[float] = 16,
    category_label_position_offsets: Optional[dict] = {},
    category_label_match_color: Optional[bool] = True,
    figsize: Optional[Iterable] = [10, 10],
    figfile: Optional[str] = None,
    show: Optional[bool] = False,
) -> Optional[mpl.axes.Axes]:
    """ """
    # process input data


circle = pycircos.Gcircle()

arc = pycircos.Garc(
    arc_id="mAbs",
    size=cc25.shape[0],
    raxis_range=(50, 100),
    linewidth=0,
    facecolor="#FFFFFF",
    label_visible=False,
)

circle.add_garc(arc)
circle.set_garcs(start=0, end=320)

for i, agbc in enumerate(["WA1", "alpha", "beta", "gamma", "kappa"]):
    vals = [1 if v else 0 for v in cc25[f"is_{agbc}"]]
    positions = list(range(len(vals)))
    widths = [0.5] * len(vals)
    raxis = [300 + 75 * i, 350 + 75 * i]
    facecolor = cdict[agbc]
    circle.barplot(
        "mAbs",
        data=vals,
        positions=positions,
        width=widths,
        rlim=[0, 1],
        raxis_range=raxis,
        facecolor=facecolor,
        edgecolor="#ffffff",
        spine=False,
    )


ax = plt.gca()
plt.text(
    0.5,
    0.515,
    "CC25",
    fontsize=36,
    weight="heavy",
    transform=ax.transAxes,
    ha="center",
    va="center",
)
plt.text(
    0.5,
    0.46,
    f"n={cc25.shape[0]}",
    fontsize=20,
    transform=ax.transAxes,
    ha="center",
    va="center",
)

plt.text(
    0.485, 0.66, "WA1", fontsize=16, transform=ax.transAxes, ha="right", va="center"
)
plt.text(
    0.485, 0.7, "alpha", fontsize=16, transform=ax.transAxes, ha="right", va="center"
)
plt.text(
    0.485, 0.74, "beta", fontsize=16, transform=ax.transAxes, ha="right", va="center"
)
plt.text(
    0.485, 0.775, "gamma", fontsize=16, transform=ax.transAxes, ha="right", va="center"
)
plt.text(
    0.485, 0.8125, "kappa", fontsize=16, transform=ax.transAxes, ha="right", va="center"
)

circle.save("./figures/CC25_agbc-ringplot_pairs-only_HSA-9999", format="pdf")
