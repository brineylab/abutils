import pytest

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import seaborn as sns

from abutils.cl import (
    get_cmap,
    truncate_colormap,
    hex_to_rgb,
    rgb_to_hex,
    hls,
    husl,
    palettes,
    show_palettes,
)


# ---------------------------
#        colormaps
# ---------------------------


def test_get_cmap_with_colormap():
    cmap = plt.get_cmap("viridis")
    assert get_cmap(cmap) == cmap


def test_get_cmap_with_colormap_name():
    cmap = plt.get_cmap("viridis")
    assert get_cmap("viridis") == cmap


def test_get_cmap_with_hex_code():
    cmap = LinearSegmentedColormap.from_list("", ["#FFFFFF", "#FF0000"])
    assert get_cmap("#FF0000") == cmap


def test_get_cmap_with_matplotlib_color():
    cmap = LinearSegmentedColormap.from_list("", ["#FFFFFF", "#0000FF"])
    assert get_cmap("blue") == cmap


def test_get_cmap_with_rgb_tuple():
    cmap = LinearSegmentedColormap.from_list("", ["#FFFFFF", "#00FF00"])
    assert get_cmap((0, 1, 0)) == cmap


def test_get_cmap_with_single_color():
    cmap = LinearSegmentedColormap.from_list("", ["#FFFFFF", "#FF0000"])
    assert get_cmap("#FF0000") == cmap


def test_get_cmap_with_dark():
    cmap = LinearSegmentedColormap.from_list("", ["#FF0000", "#000000"])
    assert get_cmap("#FF0000", dark=True) == cmap


def test_get_cmap_with_zero_color():
    cmap = LinearSegmentedColormap.from_list("", ["#C7C7C7", "#FF0000"])
    assert get_cmap("#FF0000", zero_color="#C7C7C7") == cmap


def test_get_cmap_with_minval_and_maxval():
    cmap = LinearSegmentedColormap.from_list("", ["#FFFFFF", "#FF0000"])
    assert get_cmap("#FF0000", minval=0.25, maxval=0.75) == cmap


def test_get_cmap_with_name():
    cmap = LinearSegmentedColormap.from_list("my_cmap", ["#FFFFFF", "#FF0000"])
    assert get_cmap("#FF0000", name="my_cmap") == cmap


def test_get_cmap_with_invalid_colormap():
    with pytest.raises(ValueError):
        get_cmap("invalid_colormap")


def test_truncate_colormap_with_full_range():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap)
    assert truncated_cmap.name == "viridis-trunc-0-1"
    assert truncated_cmap.colors == cmap.colors


def test_truncate_colormap_with_partial_range():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap, minval=0.25, maxval=0.75)
    assert truncated_cmap.name == "viridis-trunc-0.25-0.75"
    assert len(truncated_cmap.colors) == 256
    assert np.allclose(truncated_cmap.colors[0], cmap(0.25))
    assert np.allclose(truncated_cmap.colors[-1], cmap(0.75))


def test_truncate_colormap_with_reversed_range():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap, minval=0.75, maxval=0.25)
    assert truncated_cmap.name == "viridis-trunc-0.75-0.25"
    assert len(truncated_cmap.colors) == 256
    assert np.allclose(truncated_cmap.colors[0], cmap(0.75))
    assert np.allclose(truncated_cmap.colors[-1], cmap(0.25))


def test_truncate_colormap_with_custom_n():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap, n=10)
    assert truncated_cmap.name == "viridis-trunc-0-1"
    assert len(truncated_cmap.colors) == 10


# ---------------------------
#        HEX and RGB
# ---------------------------


def test_hex_to_rgb():
    assert hex_to_rgb("#FFFFFF") == (255, 255, 255)
    assert hex_to_rgb("#000000") == (0, 0, 0)
    assert hex_to_rgb("#FF0000") == (255, 0, 0)
    assert hex_to_rgb("#00FF00") == (0, 255, 0)
    assert hex_to_rgb("#0000FF") == (0, 0, 255)


def test_rgb_to_hex():
    assert rgb_to_hex((255, 255, 255)) == "#FFFFFF"
    assert rgb_to_hex((0, 0, 0)) == "#000000"
    assert rgb_to_hex((255, 0, 0)) == "#FF0000"
    assert rgb_to_hex((0, 255, 0)) == "#00FF00"
    assert rgb_to_hex((0, 0, 255)) == "#0000FF"


# ---------------------------
#         palettes
# ---------------------------


def test_hls():
    assert hls(3) == sns.hls_palette(3, h=0.01, l=0.6, s=0.65)
    assert hls(5, hue=0.5, lightness=0.4, saturation=0.8) == sns.hls_palette(
        5, h=0.5, l=0.4, s=0.8
    )


def test_husl():
    assert husl(3) == sns.husl_palette(3, h=0.01, s=0.9, l=0.65)
    assert husl(5, hue=0.5, saturation=0.8, lightness=0.4) == sns.husl_palette(
        5, h=0.5, s=0.8, l=0.4
    )


def test_show_palettes():
    with plt.show._enter_session():
        show_palettes(palettes)
