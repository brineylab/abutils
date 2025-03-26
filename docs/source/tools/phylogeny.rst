.. _phylogeny:


phylogeny
===================

``abutils`` provides utilities for computing and visualizing phylogenetic trees from sequence data. These utilities include
functions for tree inference using FastTree, tree visualization with baltic, and utilities for working with phylogenetic trees
in the context of antibody lineage analysis. The phylogeny module integrates sequence clustering, multiple sequence alignment,
and tree inference into a cohesive and easy-to-use interface.

All phylogeny functions are accessible through the ``abutils.tools.phylo`` module.

|  

.. csv-table:: 
   :header: "phylogeny utility", "description"
   :widths: 10, 16

   ":ref:`abutils.tools.phylo.fasttree() <fasttree-function>`", "computes a phylogenetic tree from a multiple sequence alignment using FastTree"
   ":ref:`abutils.tools.phylo.phylogeny() <phylogeny-function>`", "creates a Phylogeny object from a list of sequences or a FASTA file"
   ":ref:`abutils.tools.phylo.Phylogeny <phylogeny-class>`", "a class for representing and visualizing phylogenetic trees"

examples
---------

**computing a tree with FastTree**

Compute a phylogenetic tree from a multiple sequence alignment:

.. code-block:: python

    import abutils
    
    # Compute a tree from a FASTA alignment file
    tree_file = abutils.tools.phylo.fasttree(
        "my_alignment.fasta",
        tree_file="my_tree.newick",
        is_aa=False
    )
    
    # Or compute a tree from an alignment string and get the Newick tree as a string
    alignment_string = """>seq1
    ACGTACGTACGT
    >seq2
    ACGTACGTACGA
    >seq3
    ACATACGTACGA
    """
    tree_string = abutils.tools.phylo.fasttree(alignment_string)

|

**creating and visualizing a phylogeny**

Create a Phylogeny object and visualize the tree:

.. code-block:: python

    import abutils
    import matplotlib.pyplot as plt
    
    # Create a Phylogeny object from a list of Sequence objects
    sequences = [...]  # List of abutils.Sequence objects
    phylo = abutils.tools.phylo.phylogeny(
        sequences,
        name="my_lineage",
        cluster=True,
        clustering_threshold=0.97
    )
    
    # Plot the phylogenetic tree
    fig = plt.figure(figsize=(10, 8))
    ax = phylo.plot(
        size_multiplier=15,
        color="steelblue",
        linewidth=1.5,
        marker="o",
        marker_edgewidth=1,
        marker_edgecolor="black"
    )
    plt.tight_layout()
    plt.show()

|

**customizing phylogenetic tree visualization**

Customize the tree visualization with different marker sizes, colors, and layouts:

.. code-block:: python

    import abutils
    import matplotlib.pyplot as plt
    
    # Create a Phylogeny object
    phylo = abutils.tools.phylo.phylogeny("sequences.fasta")
    
    # Create a color mapping for specific sequences
    color_dict = {
        "seq1": "red",
        "seq2": "blue",
        "seq3": "green"
    }
    
    # Create a size mapping for specific sequences
    size_dict = {
        "seq1": 3,
        "seq2": 2,
        "seq3": 1
    }
    
    # Plot a radial tree with custom colors and sizes
    fig = plt.figure(figsize=(10, 10))
    ax = phylo.plot(
        color=color_dict,
        size=size_dict,
        radial=True,
        radial_start=0.1,
        radial_fraction=0.8,
        color_branches=True,
        marker="o",
        alpha=0.8
    )
    plt.tight_layout()
    plt.show()

| 

api
-------

.. _inference:

inference
-------------------

.. _fasttree-function:

.. autofunction:: abutils.tools.phylo.fasttree

.. autofunction:: abutils.utils.phylogeny.igphyml

.. autofunction:: abutils.utils.phylogeny.lsd


.. _drawing-trees:

drawing trees
-------------------

.. _phylogeny-function:

.. autofunction:: abutils.tools.phylo.phylogeny

.. _phylogeny-class:

.. autoclass:: abutils.tools.phylo.Phylogeny
   :members: plot, cluster, tree, root, sizes, clusters