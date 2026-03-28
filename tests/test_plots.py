"""
Smoke tests for the plotting module.

These tests verify that plotting functions run without error and return
expected types. They do not verify the visual correctness of the plots.
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


def test_bar_basic_smoke(categorical_data):
    """Smoke test: bar() with x parameter runs without error."""
    from abutils.plots.bar import bar

    ax = bar(x="category", data=categorical_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    plt.close("all")


def test_bar_with_hue_smoke(categorical_data):
    """Smoke test: bar() with hue parameter runs without error."""
    from abutils.plots.bar import bar

    ax = bar(x="category", hue="group", data=categorical_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    plt.close("all")


def test_bar_horizontal_smoke(categorical_data):
    """Smoke test: bar() with horizontal orientation runs without error."""
    from abutils.plots.bar import bar

    ax = bar(y="category", data=categorical_data, orientation="horizontal")
    assert isinstance(ax, matplotlib.axes.Axes)
    plt.close("all")


# ==========================
#    SCATTER PLOT TESTS
# ==========================


def test_scatter_basic_smoke(scatter_data):
    """Smoke test: scatter() runs without error and returns Axes."""
    from abutils.plots.scatter import scatter

    ax = scatter(x="x", y="y", data=scatter_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    plt.close("all")


def test_scatter_with_hue_smoke(scatter_data):
    """Smoke test: scatter() with hue parameter runs without error."""
    from abutils.plots.scatter import scatter

    ax = scatter(x="x", y="y", hue="category", data=scatter_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    plt.close("all")


# ==========================
#    HEATMAP TESTS
# ==========================


@pytest.mark.xfail(reason="heatmap requires specific data format via PlotData")
def test_heatmap_basic_smoke(heatmap_data):
    """Smoke test: heatmap() runs without error and returns Axes."""
    from abutils.plots.heatmap import heatmap

    # heatmap expects a simple 2D array/DataFrame of values
    ax = heatmap(data=heatmap_data.values.tolist())
    assert ax is not None
    plt.close("all")


# =====================
#    KDE PLOT TESTS
# =====================


def test_kde_basic_smoke(kde_data):
    """Smoke test: kde() runs without error and returns Axes."""
    from abutils.plots.kde import kde

    ax = kde(x="values", data=kde_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    plt.close("all")


def test_kde_with_hue_smoke(kde_data):
    """Smoke test: kde() with hue parameter runs without error."""
    from abutils.plots.kde import kde

    ax = kde(x="values", hue="category", data=kde_data)
    assert isinstance(ax, matplotlib.axes.Axes)
    plt.close("all")


# ======================
#    RIDGE PLOT TESTS
# ======================


def test_ridge_basic_smoke(kde_data):
    """Smoke test: ridge() runs without error."""
    from abutils.plots.ridge import ridge

    # ridge uses 'categories' and 'values' parameters
    ax = ridge(categories="category", values="values", data=kde_data, show=False)
    # ridge may return axes or figure depending on implementation
    assert ax is not None
    plt.close("all")


# ======================
#    DONUT PLOT TESTS
# ======================


def test_donut_basic_smoke(donut_data):
    """Smoke test: donut() runs without error."""
    from abutils.plots.donut import donut

    # donut uses 'values' for categories and 'counts' for values
    ax = donut(values="category", counts="value", data=donut_data)
    assert ax is not None
    plt.close("all")
