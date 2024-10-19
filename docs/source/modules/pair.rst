.. _pair-class:

pair
===================

``Pair`` objects are used to represent one or ``Sequence`` objects that together
represent a single mAb or TCR. 

|  

instantiation
----------------

An individual ``Pair`` object can be instantiated directly from a 
list of ``Sequence`` objects that belong to the same mAb or TCR:

.. code-block:: python

    import abutils

    # create a pair from a list of sequences
    pair = abutils.Pair([heavy_sequence, light_sequence])

|  

Or, more commonly, from a larger list of ``Sequence`` objects for which the appropriate pairing
needs to be determined from the sequence names, using :ref:`assign_pairs() <assign-pairs>`:

.. code-block:: python

    import abutils

    # batch create pairs from a list of annotated sequences
    sequences = abutils.io.read_airr('path/to/sequences.tsv')
    pairs = abutils.tl.assign_pairs(sequences)

|  

api
----------------
.. autoclass:: abutils.core.pair.Pair

.. _assign-pairs:

.. autofunction:: abutils.core.pair.assign_pairs


