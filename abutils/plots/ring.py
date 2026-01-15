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
from typing import Iterable

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

from natsort import natsorted

import pycircos


def ring(
    values: str | dict | Iterable,
    categories: str | Iterable | None = None,
    data: pd.DataFrame | None = None,
    color: str | Iterable | None = None,
    palette: str | Iterable | dict | None = None,
    raxis_range: tuple | None = None,
    linewidth: float | None = 1.0,
    facecolor: str | None = "#FFFFFF",
    label_visible: bool | None = False,
    arc_start: float | None = 0,
    arc_end: float | None = 360,
    title: str | None = None,
    title_fontsize: float | None = 36,
    title_fontweight: str | None = "heavy",
    subtitle: str | None = None,
    subtitle_fontsize: float | None = 20,
    subtitle_fontweight: str | None = "normal",
    show_category_labels: bool | None = True,
    category_label_fontsize: float | None = 16,
    category_label_position_offsets: dict | None = {},
    category_label_match_color: bool | None = True,
    figsize: Iterable | None = [10, 10],
    figfile: str | None = None,
    show: bool | None = False,
) -> mpl.axes.Axes | None:
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
