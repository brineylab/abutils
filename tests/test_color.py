import pytest

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import seaborn as sns

from ..utils.color import (
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
    got_cmap = get_cmap(cmap)
    for i in np.arange(0, 1, 0.1):
        assert np.allclose(got_cmap(i), cmap(i))


def test_get_cmap_with_colormap_name():
    cmap = plt.get_cmap("viridis")
    got_cmap = get_cmap("viridis")
    for i in np.arange(0, 1, 0.1):
        assert np.allclose(got_cmap(i), cmap(i))


def test_get_cmap_with_hex_code():
    cmap = sns.light_palette("#FF0000", as_cmap=True)
    got_cmap = get_cmap("#FF0000")
    for i in np.arange(0, 1, 0.1):
        assert np.allclose(got_cmap(i), cmap(i))


def test_get_cmap_with_matplotlib_color():
    cmap = sns.light_palette("#0000FF", as_cmap=True)
    got_cmap = get_cmap("blue")
    for i in np.arange(0, 1, 0.1):
        assert np.allclose(got_cmap(i), cmap(i))


def test_get_cmap_with_rgb_tuple():
    cmap = sns.light_palette("#00FF00", as_cmap=True)
    got_cmap = get_cmap((0, 1, 0))
    for i in np.arange(0, 1, 0.1):
        assert np.allclose(got_cmap(i), cmap(i))


def test_get_cmap_with_single_color():
    cmap = sns.light_palette("#FF0000", as_cmap=True)
    got_cmap = get_cmap("#FF0000")
    for i in np.arange(0, 1, 0.1):
        assert np.allclose(got_cmap(i), cmap(i))


def test_get_cmap_with_dark():
    cmap = sns.dark_palette("#FF0000", as_cmap=True)
    got_cmap = get_cmap("#FF0000", dark=True)
    for i in np.arange(0, 1, 0.1):
        assert np.allclose(got_cmap(i), cmap(i))


def test_get_cmap_with_zero_color():
    got_cmap = get_cmap("#FF0000", zero_color=(0, 0, 0))
    assert got_cmap(0) == (0, 0, 0, 1)


def test_get_cmap_with_minval_and_maxval():
    cmap = get_cmap("#FF0000", name="my_cmap")
    got_cmap = get_cmap("#FF0000", minval=0.25, maxval=0.75, name="truncated_cmap")
    assert np.allclose(got_cmap(0), cmap(0.25))
    assert np.allclose(got_cmap.reversed()(0), cmap.reversed()(0.25), rtol=0.05)
    # assert np.allclose(got_cmap(1), cmap(0.75))


def test_get_cmap_with_name():
    cmap = get_cmap("#FF0000", name="my_cmap")
    cmap.name = "my_cmap"


@pytest.mark.xfail
def test_get_cmap_with_invalid_colormap():
    with pytest.raises(ValueError):
        get_cmap("invalid_colormap")


def test_truncate_colormap_with_full_range():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap)
    assert truncated_cmap.name == "viridis-trunc-0-1"
    assert truncated_cmap.N == 256
    assert np.allclose(truncated_cmap(0), cmap(0))
    assert np.allclose(truncated_cmap(1), cmap(1))


def test_truncate_colormap_with_partial_range():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap, minval=0.25, maxval=0.75)
    assert truncated_cmap.name == "viridis-trunc-0.25-0.75"
    assert truncated_cmap.N == 256
    assert np.allclose(truncated_cmap(0), cmap(0.25))
    assert np.allclose(truncated_cmap.reversed()(0), cmap(0.75), rtol=0.05)


def test_truncate_colormap_with_reversed_range():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap, minval=0.75, maxval=0.25)
    assert truncated_cmap.name == "viridis-trunc-0.75-0.25"
    assert truncated_cmap.N == 256
    assert np.allclose(truncated_cmap(0), cmap(0.75))
    assert np.allclose(truncated_cmap.reversed()(1), cmap(0.25), rtol=0.05)


def test_truncate_colormap_with_custom_n():
    cmap = plt.get_cmap("viridis")
    truncated_cmap = truncate_colormap(cmap, n=10)
    assert truncated_cmap.name == "viridis-trunc-0-1"
    assert truncated_cmap.N == 10


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
    assert rgb_to_hex((255, 255, 255)) == "#ffffff"
    assert rgb_to_hex((0, 0, 0)) == "#000000"
    assert rgb_to_hex((255, 0, 0)) == "#ff0000"
    assert rgb_to_hex((0, 255, 0)) == "#00ff00"
    assert rgb_to_hex((0, 0, 255)) == "#0000ff"


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


@pytest.mark.skip(reason="Not sure how to test plotting functions")
def test_show_palettes():
    with plt.show._enter_session():
        show_palettes(palettes)
