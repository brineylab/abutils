.. _kde-plot:

kde plot
========

``abutils`` provides functionality for creating kernel density estimate (KDE) plots. These plots are useful
for visualizing the distribution of continuous data, either in one or two dimensions. The KDE function combines
the flexibility of scatter plots with density visualization, allowing for the integration of categorical or
continuous hue variables.

|  

.. csv-table:: 
   :header: "plotting function", "description"
   :widths: 10, 16

   ":ref:`abutils.plots.kde() <kde-function>`", "creates KDE plots with optional scatter overlay and hue coloring"

examples
---------

**basic KDE plot**

A simple one-dimensional KDE plot shows the distribution of a single variable:

.. code-block:: python

    import abutils
    import numpy as np

    # create sample data
    x = np.random.normal(0, 1, 1000)
    
    # create a 1D KDE plot
    ax = abutils.plots.kde(
        x=x,
        fill=True,
        xlabel='Value',
        ylabel='Density',
        figsize=[8, 5]
    )

|

**two-dimensional KDE with scatter points**

Adding a y-parameter creates a two-dimensional KDE plot, and points can be overlaid:

.. code-block:: python

    import abutils
    import numpy as np
    import pandas as pd

    # create correlated data
    x = np.random.normal(0, 1, 500)
    y = x + np.random.normal(0, 0.5, 500)
    
    # create a 2D KDE plot with scatter points
    ax = abutils.plots.kde(
        x=x,
        y=y,
        show_scatter=True,
        scatter_size=10,
        scatter_alpha=0.5,
        fill=True,
        kde_fill_alpha=0.5,
        xlabel='Variable X',
        ylabel='Variable Y',
        figsize=[8, 8]
    )

|

**KDE with categorical hue**

You can color the KDE and scatter points by a categorical variable:

.. code-block:: python

    import abutils
    import numpy as np
    import pandas as pd

    # create data with categories
    n_points = 600
    data = pd.DataFrame({
        'x': np.concatenate([np.random.normal(-1, 0.5, n_points//3), 
                            np.random.normal(0, 0.5, n_points//3),
                            np.random.normal(1, 0.5, n_points//3)]),
        'y': np.concatenate([np.random.normal(1, 0.5, n_points//3), 
                            np.random.normal(0, 0.5, n_points//3),
                            np.random.normal(-1, 0.5, n_points//3)]),
        'category': np.repeat(['A', 'B', 'C'], n_points//3)
    })
    
    # create a KDE plot with categorical hue
    ax = abutils.plots.kde(
        x='x',
        y='y',
        hue='category',
        data=data,
        palette={'A': 'red', 'B': 'blue', 'C': 'green'},
        fill=True,
        kde_fill_alpha=0.3,
        scatter_alpha=0.5,
        scatter_size=15,
        xlabel='X Value',
        ylabel='Y Value',
        figsize=[10, 8]
    )

|

**KDE with continuous hue and highlighted points**

KDE plots support continuous hue values and allow highlighting specific points:

.. code-block:: python

    import abutils
    import numpy as np
    import pandas as pd

    # create data with continuous hue values
    n_points = 500
    x = np.random.normal(0, 1, n_points)
    y = x + np.random.normal(0, 0.5, n_points)
    hue_val = np.abs(x + y)  # continuous hue based on distance from origin
    
    # identify points to highlight
    highlight_indices = np.where(hue_val > 2.5)[0]
    
    # create a KDE plot with continuous hue and highlighted points
    ax = abutils.plots.kde(
        x=x,
        y=y,
        hue=hue_val,
        highlight_index=highlight_indices,
        highlight_marker='*',
        highlight_size=100,
        highlight_color='black',
        highlight_name='Extreme Values',
        cmap='viridis',
        scatter_size=20,
        scatter_alpha=0.7,
        equal_axes=True,
        xlabel='X Value',
        ylabel='Y Value',
        figsize=[10, 8]
    )

| 

api
-------

.. _kde-function:

.. autofunction:: abutils.plots.kde 