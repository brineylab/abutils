.. _scatter-plot:

scatter
===========

``abutils`` provides functionality for creating highly customizable scatter plots. The scatter plot function
supports a wide range of features including categorical and continuous color scaling, point highlighting,
and legend customization. It is designed to work seamlessly with data in various formats, including
pandas DataFrames and abutils Sequence objects.

|  

.. csv-table:: 
   :header: "plotting function", "description"
   :widths: 10, 16

   ":ref:`abutils.plots.scatter() <scatter-function>`", "creates scatter plots with extensive customization options"

examples
---------

**basic scatter plot**

Create a simple scatter plot from x and y coordinates:

.. code-block:: python

    import abutils
    import numpy as np

    # create sample data
    x = np.random.normal(0, 1, 100)
    y = x + np.random.normal(0, 0.5, 100)
    
    # create a basic scatter plot
    ax = abutils.plots.scatter(
        x=x,
        y=y,
        color='steelblue',
        size=50,
        alpha=0.7,
        xlabel='X Variable',
        ylabel='Y Variable',
        figsize=[8, 8]
    )

|

**scatter plot with categorical hue**

Color points based on a categorical variable:

.. code-block:: python

    import abutils
    import numpy as np
    import pandas as pd

    # create data with categories
    n_points = 150
    data = pd.DataFrame({
        'x': np.concatenate([np.random.normal(-1, 0.5, n_points//3), 
                            np.random.normal(0, 0.5, n_points//3),
                            np.random.normal(1, 0.5, n_points//3)]),
        'y': np.concatenate([np.random.normal(1, 0.5, n_points//3), 
                            np.random.normal(0, 0.5, n_points//3),
                            np.random.normal(-1, 0.5, n_points//3)]),
        'group': np.repeat(['Group A', 'Group B', 'Group C'], n_points//3)
    })
    
    # create a scatter plot with categorical hue
    ax = abutils.plots.scatter(
        x='x',
        y='y',
        hue='group',
        data=data,
        palette={'Group A': 'crimson', 'Group B': 'navy', 'Group C': 'forestgreen'},
        size=80,
        alpha=0.7,
        xlabel='X Coordinate',
        ylabel='Y Coordinate',
        equal_axes=True,
        figsize=[10, 8]
    )

|

**scatter plot with continuous hue and color bar**

Color points based on a continuous variable with a custom colormap:

.. code-block:: python

    import abutils
    import numpy as np
    import pandas as pd

    # create data with a continuous variable
    x = np.random.uniform(-3, 3, 200)
    y = np.random.uniform(-3, 3, 200)
    
    # calculate distance from origin as the hue value
    distance = np.sqrt(x**2 + y**2)
    
    # create a scatter plot with continuous hue
    ax = abutils.plots.scatter(
        x=x,
        y=y,
        hue=distance,
        cmap='plasma',  # custom colormap
        size=70,
        alpha=0.8,
        xlabel='X Coordinate',
        ylabel='Y Coordinate',
        cbar_title='Distance from Origin',
        cbar_title_fontsize=14,
        equal_axes=True,
        figsize=[10, 8]
    )

|

**scatter plot with highlighted points**

Highlight specific points of interest in the scatter plot:

.. code-block:: python

    import abutils
    import numpy as np
    import pandas as pd

    # create random data
    np.random.seed(42)
    n_points = 100
    x = np.random.normal(0, 1, n_points)
    y = np.random.normal(0, 1, n_points)
    
    # identify points to highlight (those far from the origin)
    distance = np.sqrt(x**2 + y**2)
    highlight_indices = np.where(distance > 2.0)[0]
    
    # create a scatter plot with highlighted points
    ax = abutils.plots.scatter(
        x=x,
        y=y,
        highlight_index=highlight_indices,
        highlight_marker='*',
        highlight_size=200,
        highlight_color='red',
        highlight_name='Outliers',
        color='darkblue',
        size=50,
        alpha=0.6,
        xlabel='X Value',
        ylabel='Y Value',
        equal_axes=True,
        figsize=[10, 8]
    )

|

**scatter plot with customized legend**

Create a scatter plot with a custom-positioned legend:

.. code-block:: python

    import abutils
    import numpy as np
    import pandas as pd

    # create data with categorical and continuous variables
    n_points = 200
    categories = ['Low', 'Medium', 'High']
    data = pd.DataFrame({
        'x': np.random.normal(0, 1, n_points),
        'y': np.random.normal(0, 1, n_points),
        'category': np.random.choice(categories, n_points),
        'value': np.random.uniform(0, 10, n_points)
    })
    
    # create a scatter plot with custom legend
    ax = abutils.plots.scatter(
        x='x',
        y='y',
        hue='category',
        data=data,
        size=80,
        alpha=0.7,
        palette={'Low': 'blue', 'Medium': 'green', 'High': 'red'},
        legend_on_data=True,  # place legend on the data
        legend_marker_alpha=1.0,  # make legend markers fully opaque
        legend_fontsize=14,
        legend_fontweight='bold',
        xlabel='X Variable',
        ylabel='Y Variable',
        figsize=[12, 10]
    )

| 

api
-------

.. _scatter-function:

.. autofunction:: abutils.plots.scatter 