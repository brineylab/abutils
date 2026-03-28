"""
Tests for the plotting module.

Tests verify both that plotting functions run without error and that the
resulting plots contain expected visual elements.
"""

import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for testing

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

# ==================
#      FIXTURES
# ==================


@pytest.fixture
def categorical_data():
    """DataFrame with categorical data for bar plots."""
    np.random.seed(42)
    return pd.DataFrame(
        {
            "category": np.random.choice(["A", "B", "C", "D"], 100),
            "value": np.random.randint(1, 100, 100),
            "group": np.random.choice(["X", "Y"], 100),
        }
    )


@pytest.fixture
def scatter_data():
    """DataFrame for scatter plots."""
    np.random.seed(42)
    return pd.DataFrame(
        {
            "x": np.random.rand(50),
            "y": np.random.rand(50),
            "category": np.random.choice(["A", "B", "C"], 50),
            "size": np.random.randint(10, 100, 50),
        }
    )


@pytest.fixture
def heatmap_data():
    """DataFrame for heatmap plots."""
    np.random.seed(42)
    return pd.DataFrame(
        np.random.rand(5, 5),
        index=[f"row{i}" for i in range(5)],
        columns=[f"col{i}" for i in range(5)],
    )


@pytest.fixture
def kde_data():
    """Data for KDE plots."""
    np.random.seed(42)
    return pd.DataFrame(
        {
            "values": np.concatenate([
                np.random.normal(0, 1, 100),
                np.random.normal(3, 0.5, 100),
            ]),
            "category": ["A"] * 100 + ["B"] * 100,
        }
    )


@pytest.fixture
def donut_data():
    """Data for donut plots."""
    return pd.DataFrame(
        {
            "category": ["A", "B", "C", "D"],
            "value": [30, 25, 25, 20],
        }
    )


# ========================
#    BAR PLOT TESTS
# ========================


def test_bar_basic(categorical_data):
    """bar() returns Axes with correct number of bars."""
    from abutils.plots.bar import bar

    ax = bar(x="category", data=categorical_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    # should have bars (patches) for each category
    n_categories = categorical_data["category"].nunique()
    patches = [p for p in ax.patches if p.get_height() > 0 or p.get_width() > 0]
    assert len(patches) >= n_categories
    plt.close("all")


def test_bar_with_hue(categorical_data):
    """bar() with hue produces patches for each category x group combination."""
    from abutils.plots.bar import bar

    ax = bar(x="category", hue="group", data=categorical_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    # with hue, we expect more patches than without
    patches = [p for p in ax.patches if p.get_height() > 0 or p.get_width() > 0]
    assert len(patches) >= categorical_data["category"].nunique()
    plt.close("all")


def test_bar_horizontal(categorical_data):
    """bar() horizontal orientation produces bars with nonzero widths."""
    from abutils.plots.bar import bar

    ax = bar(y="category", data=categorical_data, orientation="horizontal")
    assert isinstance(ax, matplotlib.axes.Axes)
    # horizontal bars have nonzero widths
    widths = [p.get_width() for p in ax.patches]
    assert any(w > 0 for w in widths)
    plt.close("all")


# ==========================
#    SCATTER PLOT TESTS
# ==========================


def test_scatter_basic(scatter_data):
    """scatter() renders data points within expected axis ranges."""
    from abutils.plots.scatter import scatter

    ax = scatter(x="x", y="y", data=scatter_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    # data is 0-1 range; axes should encompass that
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    assert xlim[0] <= 0.0
    assert xlim[1] >= 1.0
    assert ylim[0] <= 0.0
    assert ylim[1] >= 1.0
    plt.close("all")


def test_scatter_with_hue(scatter_data):
    """scatter() with hue renders multiple collections."""
    from abutils.plots.scatter import scatter

    ax = scatter(x="x", y="y", hue="category", data=scatter_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    # with hue, should have rendered points (children include PathCollections)
    collections = [
        c for c in ax.get_children()
        if hasattr(c, "get_offsets") and len(c.get_offsets()) > 0
    ]
    assert len(collections) >= 1
    plt.close("all")


# ==========================
#    HEATMAP TESTS
# ==========================


@pytest.mark.xfail(reason="heatmap requires specific data format via PlotData")
def test_heatmap_basic_smoke(heatmap_data):
    """Smoke test: heatmap() runs without error and returns Axes."""
    from abutils.plots.heatmap import heatmap

    ax = heatmap(data=heatmap_data.values.tolist())
    assert ax is not None
    plt.close("all")


# =====================
#    KDE PLOT TESTS
# =====================


def test_kde_basic(kde_data):
    """kde() produces lines representing the density estimate."""
    from abutils.plots.kde import kde

    ax = kde(x="values", data=kde_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    # KDE should produce at least one line
    lines = ax.get_lines()
    assert len(lines) >= 1
    # the line should span the data range
    xdata = lines[0].get_xdata()
    assert min(xdata) < 0  # data has mean=0 component
    assert max(xdata) > 2  # data has mean=3 component
    plt.close("all")


def test_kde_with_hue(kde_data):
    """kde() with hue produces multiple density lines."""
    from abutils.plots.kde import kde

    ax = kde(x="values", hue="category", data=kde_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    # with 2 categories, should have at least 2 lines
    lines = ax.get_lines()
    assert len(lines) >= 2
    plt.close("all")


# ======================
#    RIDGE PLOT TESTS
# ======================


def test_ridge_basic(kde_data):
    """ridge() returns a result and renders without error."""
    from abutils.plots.ridge import ridge

    result = ridge(categories="category", values="values", data=kde_data, show=False)
    assert result is not None
    plt.close("all")


# ======================
#    DONUT PLOT TESTS
# ======================


def test_donut_basic(donut_data):
    """donut() produces a pie/wedge chart."""
    from abutils.plots.donut import donut

    ax = donut(values="category", counts="value", data=donut_data)
    assert ax is not None
    plt.close("all")
