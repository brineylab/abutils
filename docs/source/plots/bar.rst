.. _bar-plot:

bar plot
========

``abutils`` provides functionality for creating bar plots of categorical data. The bar plot function
is flexible and can handle both raw counts and frequency distributions, with options for highlighting
specific categories, customizing colors, and adding legends.

|  

.. csv-table:: 
   :header: "plotting function", "description"
   :widths: 10, 16

   ":ref:`abutils.plots.bar() <bar-function>`", "creates bar plots for categorical data with optional grouping"

examples
---------

**basic bar plot**

Bar plots can be created using either raw values or named columns in a DataFrame. If only x-values are 
provided, the function counts occurrences automatically.

.. code-block:: python

    import abutils

    # create a simple bar plot from a list of categories
    categories = ['A', 'B', 'A', 'C', 'B', 'A', 'D', 'A', 'B', 'C']
    ax = abutils.plots.bar(x=categories, figsize=[8, 5])

|

**stacked bar plot with hue**

Adding a hue parameter creates a stacked bar plot where each segment represents a hue category.

.. code-block:: python

    import abutils
    import pandas as pd

    # create a DataFrame with categories and subcategories
    data = pd.DataFrame({
        'gene': ['IGHV1-2', 'IGHV1-2', 'IGHV1-18', 'IGHV3-23', 'IGHV3-23', 'IGHV4-34'],
        'isotype': ['IgG', 'IgM', 'IgG', 'IgA', 'IgG', 'IgM']
    })

    # create a stacked bar plot
    ax = abutils.plots.bar(
        x='gene', 
        hue='isotype', 
        data=data,
        palette={'IgG': 'darkred', 'IgM': 'navy', 'IgA': 'forestgreen'},
        figsize=[10, 5]
    )

|

**normalized frequency bar plot**

Bar plots can show normalized frequencies instead of raw counts, and can be displayed horizontally.

.. code-block:: python

    import abutils
    import pandas as pd

    # create a DataFrame with counts
    data = pd.DataFrame({
        'category': ['A', 'B', 'C', 'D', 'E'],
        'count': [23, 45, 12, 67, 8]
    })

    # create a horizontal bar plot with normalized frequencies
    ax = abutils.plots.bar(
        x='category',
        y='count',
        data=data,
        normalize=True,
        orientation='horizontal',
        color='steelblue',
        xlabel='Frequency (%)',
        ylabel='Category',
        figsize=[6, 8]
    )

|

**highlighting specific categories**

You can highlight specific categories to draw attention to them.

.. code-block:: python

    import abutils
    import pandas as pd

    # create a DataFrame with categories and values
    data = pd.DataFrame({
        'gene': ['IGHV1-2', 'IGHV1-18', 'IGHV3-23', 'IGHV3-30', 'IGHV4-34', 'IGHV5-51'],
        'frequency': [0.15, 0.08, 0.25, 0.12, 0.30, 0.10]
    })

    # create a bar plot with highlighted categories
    ax = abutils.plots.bar(
        x='gene',
        y='frequency',
        data=data,
        highlight=['IGHV3-23', 'IGHV4-34'],
        highlight_color='darkred',
        alt_color='lightgray',
        ylabel='Frequency',
        figsize=[10, 5]
    )

| 

api
-------

.. _bar-function:

.. autofunction:: abutils.plots.bar 