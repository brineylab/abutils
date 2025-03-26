.. _colors:

colors
======

``abutils`` provides a set of utilities for working with colors in visualizations. These utilities include
functions for creating and manipulating colormaps, generating color palettes, and converting between
different color representations. The module also includes a collection of predefined color palettes
designed for use in scientific visualizations.

All color functions are accessible through the ``abutils.cl`` module, which is the recommended way to 
use these utilities in your code.

|  

.. csv-table:: 
   :header: "color utility", "description"
   :widths: 10, 16

   ":ref:`abutils.cl.get_cmap() <get-cmap-function>`", "creates or retrieves a Matplotlib colormap"
   ":ref:`abutils.cl.monochrome_palette() <monochrome-palette-function>`", "creates a monochromatic color palette"
   ":ref:`abutils.cl.truncate_colormap() <truncate-colormap-function>`", "creates a subset of an existing colormap"
   ":ref:`abutils.cl.hex_to_rgb() <hex-to-rgb-function>`", "converts a hex color code to RGB values"
   ":ref:`abutils.cl.rgb_to_hex() <rgb-to-hex-function>`", "converts RGB values to a hex color code"
   ":ref:`abutils.cl.hls() <hls-function>`", "creates a color palette using the HLS color space"
   ":ref:`abutils.cl.husl() <husl-function>`", "creates a color palette using the HUSL color space"
   ":ref:`abutils.cl.show_palettes() <show-palettes-function>`", "displays the predefined color palettes"
   ":ref:`abutils.cl.palettes <palettes-constant>`", "dictionary of predefined color palettes"
   ":ref:`abutils.cl.true_false_palette <true-false-palette-constant>`", "a predefined palette for boolean values"

examples
---------

**using predefined color palettes**

``abutils`` includes several predefined color palettes that are designed for scientific visualizations:

.. code-block:: python

    import abutils
    import matplotlib.pyplot as plt
    
    # get a predefined palette
    palette = abutils.cl.palettes['vibrant']
    
    # use the palette in a plot
    for i, color in enumerate(palette):
        plt.plot([0, 1], [i, i], color=color, linewidth=10)
        
    plt.ylim(-0.5, len(palette) - 0.5)
    plt.title('Vibrant Palette')
    plt.show()
    
    # view all available palettes
    abutils.cl.show_palettes()

|

**creating custom colormaps**

Create custom colormaps from a single color or modify existing colormaps:

.. code-block:: python

    import abutils
    import numpy as np
    import matplotlib.pyplot as plt
    
    # create a colormap from a single color
    cmap1 = abutils.cl.get_cmap('steelblue')
    
    # create a colormap with a special color for zero values
    cmap2 = abutils.cl.get_cmap('viridis', zero_color='lightgray')
    
    # truncate an existing colormap
    cmap3 = abutils.cl.truncate_colormap(plt.cm.plasma, minval=0.2, maxval=0.8)
    
    # display the colormaps
    fig, axes = plt.subplots(3, 1, figsize=(8, 4))
    
    for ax, cmap, title in zip(axes, [cmap1, cmap2, cmap3], 
                              ['Monochrome', 'Zero-colored', 'Truncated']):
        gradient = np.linspace(0, 1, 256)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(gradient, aspect='auto', cmap=cmap)
        ax.set_title(title)
        ax.set_axis_off()
    
    plt.tight_layout()
    plt.show()

|

**generating color palettes**

Generate color palettes with different color generation methods:

.. code-block:: python

    import abutils
    import matplotlib.pyplot as plt
    
    # create a monochromatic palette
    mono_palette = abutils.cl.monochrome_palette('darkred', n_colors=7)
    
    # create an HLS palette
    hls_palette = abutils.cl.hls(7, hue=0.6, lightness=0.6, saturation=0.8)
    
    # create a HUSL palette
    husl_palette = abutils.cl.husl(7, hue=0.1, saturation=0.9, lightness=0.65)
    
    # display the palettes
    fig, axes = plt.subplots(3, 1, figsize=(8, 3))
    
    for ax, palette, title in zip(axes, [mono_palette, hls_palette, husl_palette], 
                                ['Monochrome', 'HLS', 'HUSL']):
        for i, color in enumerate(palette):
            ax.plot([i, i+0.9], [0, 0], color=color, linewidth=20)
        ax.set_xlim(-0.1, len(palette))
        ax.set_ylim(-0.5, 0.5)
        ax.set_title(title)
        ax.set_axis_off()
    
    plt.tight_layout()
    plt.show()

|

**color conversion**

Convert between RGB and hex color representations:

.. code-block:: python

    import abutils
    
    # convert hex to RGB
    rgb = abutils.cl.hex_to_rgb('#3498db')
    print(f"Hex #3498db as RGB: {rgb}")
    
    # convert RGB to hex
    hex_code = abutils.cl.rgb_to_hex((52, 152, 219))
    print(f"RGB (52, 152, 219) as hex: {hex_code}")
    
    # convert normalized RGB (0-1) to hex
    hex_code_norm = abutils.cl.rgb_to_hex((0.2, 0.6, 0.85))
    print(f"Normalized RGB (0.2, 0.6, 0.85) as hex: {hex_code_norm}")

| 

api
-------

.. _get-cmap-function:

.. autofunction:: abutils.cl.get_cmap

.. _monochrome-palette-function:

.. autofunction:: abutils.cl.monochrome_palette

.. _truncate-colormap-function:

.. autofunction:: abutils.cl.truncate_colormap

.. _hex-to-rgb-function:

.. autofunction:: abutils.cl.hex_to_rgb

.. _rgb-to-hex-function:

.. autofunction:: abutils.cl.rgb_to_hex

.. _hls-function:

.. autofunction:: abutils.cl.hls

.. _husl-function:

.. autofunction:: abutils.cl.husl

.. _show-palettes-function:

.. autofunction:: abutils.cl.show_palettes

.. _palettes-constant:

.. autodata:: abutils.cl.palettes
   :annotation: = dictionary of predefined color palettes

.. _true-false-palette-constant:

.. autodata:: abutils.cl.true_false_palette
   :annotation: = {True: "#e41a1c", False: "#d1d1d1"}
