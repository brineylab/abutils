.. _donut-plot:

donut plot
==========

``abutils`` provides a function for creating donut plots, which are useful for visualizing the distribution
of categorical data. Donut plots are particularly effective for showing the relative proportions of different
categories within a dataset, with arc widths proportional to category size. The donut plot function supports
a variety of customization options including hue coloring, title customization, and singleton grouping.

|  

.. csv-table:: 
   :header: "plotting function", "description"
   :widths: 10, 16

   ":ref:`abutils.plots.donut() <donut-function>`", "creates donut plots for categorical data with optional hue coloring"

examples
---------

**basic donut plot**

Create a simple donut plot from a list of categories:

.. code-block:: python

    import abutils

    # create a list of categories
    categories = ['A', 'A', 'A', 'B', 'B', 'C', 'C', 'C', 'C', 'D', 'E', 'E']
    
    # create a basic donut plot
    ax = abutils.plots.donut(
        values=categories,
        title='Category Distribution',
        figsize=[8, 8]
    )

|

**donut plot with custom colors and sorting**

You can customize colors and determine how segments are sorted:

.. code-block:: python

    import abutils
    import pandas as pd

    # create a DataFrame with categories and counts
    data = pd.DataFrame({
        'gene': ['IGHV1-2', 'IGHV1-18', 'IGHV3-23', 'IGHV3-30', 'IGHV4-34', 'IGHV5-51'],
        'count': [150, 80, 250, 120, 300, 100]
    })
    
    # create a donut plot with custom colors, sorted by count
    ax = abutils.plots.donut(
        values='gene',
        counts='count',
        data=data,
        color='steelblue',  # base color for monochromatic palette
        sort_by='count',
        sort_descending=True,
        title='V Gene Usage',
        subtitle='n=1000',
        figsize=[10, 10]
    )

|

**donut plot with categorical hue**

Color the segments based on a categorical hue variable:

.. code-block:: python

    import abutils
    import pandas as pd

    # create a DataFrame with categories, counts, and a hue variable
    data = pd.DataFrame({
        'lineage': ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8'],
        'count': [45, 33, 28, 20, 15, 12, 8, 5],
        'isotype': ['IgG', 'IgM', 'IgG', 'IgA', 'IgG', 'IgM', 'IgA', 'IgG']
    })
    
    # create a donut plot with categorical hue
    ax = abutils.plots.donut(
        values='lineage',
        counts='count',
        hue='isotype',
        data=data,
        palette={'IgG': 'darkred', 'IgM': 'navy', 'IgA': 'forestgreen'},
        sort_by='count',
        sort_descending=True,
        title='Lineage Distribution',
        subtitle='by Isotype',
        width=0.6,  # width of the donut
        figsize=[10, 10]
    )

|

**donut plot with grouped singletons and continuous hue**

Group small segments as singletons and color by a continuous variable:

.. code-block:: python

    import abutils
    import pandas as pd
    import numpy as np

    # create a larger dataset with many small categories
    np.random.seed(42)
    n_categories = 50
    categories = [f'Cat{i}' for i in range(1, n_categories+1)]
    
    # create a distribution with many small values (power law)
    alpha = 1.5
    counts = np.random.power(alpha, n_categories) * 1000
    counts = counts.astype(int) + 1  # ensure no zeros
    
    # create some continuous values for hue
    hue_values = np.random.uniform(0, 100, n_categories)
    
    data = pd.DataFrame({
        'category': categories,
        'count': counts,
        'metric': hue_values
    })
    
    # create a donut plot with grouped singletons and continuous hue
    ax = abutils.plots.donut(
        values='category',
        counts='count',
        hue='metric',
        data=data,
        group_singletons=True,  # group categories with count=1
        singleton_color='lightgray',
        cmap='viridis',  # colormap for continuous hue
        sort_by='count',
        sort_descending=True,
        force_continuous_hue=True,
        title=f'Categories (n={data["count"].sum()})',
        subtitle=f'Showing {n_categories} categories',
        width=0.55,
        figsize=[12, 12]
    )

| 

api
-------

.. _donut-function:

.. autofunction:: abutils.plots.donut 