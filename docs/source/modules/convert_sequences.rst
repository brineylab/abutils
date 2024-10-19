

convert sequence data
==============================


``abutils`` provides functions for converting annotated sequence data between ``Sequence`` and 
``Pair`` objects and Polars_ and Pandas_ DataFrames. All annotations are assumed to be in 
AIRR-C_ format.

|  

.. table::
   :align: left
   :widths: 10, 14, 24
   :width: 100%

   +----------+-------------------------------------------------+----------------------------------------+
   | format   | function                                        | notes                                  |
   +==========+=================================================+========================================+
   | Pandas   | :ref:`abutils.io.to_pandas() <to-pandas>`       | returns a Pandas DataFrame             |
   +          +-------------------------------------------------+----------------------------------------+
   |          | :ref:`abutils.io.from_pandas() <from-pandas>`   | returns a list of ``Sequence`` or      |
   |          |                                                 | ``Pair`` objects                       |
   +----------+-------------------------------------------------+----------------------------------------+
   | Polars   | :ref:`abutils.io.to_polars() <to-polars>`       | returns a Polars DataFrame             |
   +          +-------------------------------------------------+----------------------------------------+
   |          | :ref:`abutils.io.from_polars() <from-polars>`   | returns a list of ``Sequence`` or      |
   |          |                                                 | ``Pair`` objects                       |
   +----------+-------------------------------------------------+----------------------------------------+

|  
Pandas
------------------

``abutils`` can convert lists of ``Sequence`` or ``Pair`` objects to and from Pandas DataFrames:

.. warning::

    While ``abutils`` can convert between ``Sequence`` and ``Pair`` objects and Pandas DataFrames, the input must contains
    only one type of object. For example, you cannot mix ``Sequence`` and ``Pair`` objects in the 
    same list.

.. code-block:: python

    # convert list of sequences to Pandas DataFrame
    sequences_df = abutils.io.to_pandas(
        sequences, 
    )

    # convert Pandas DataFrame back to list of sequences
    sequences = abutils.io.from_pandas(sequences_df)

    # convert list of pairs to Pandas DataFrame
    pairs_df = abutils.io.to_pandas(
        pairs, 
    )

    # convert Pandas DataFrame back to list of pairs
    pairs = abutils.io.from_pandas(pairs_df)

|  

Polars
------------------------------
``abutils`` can convert lists of ``Sequence`` or ``Pair`` objects to and from Pandas DataFrames:

.. code-block:: python

    # convert list of sequences to Polars DataFrame
    sequences_df = abutils.io.to_polars(
        sequences, 
    )

    # convert Polars DataFrame back to list of sequences
    sequences = abutils.io.from_polars(sequences_df)

    # convert list of pairs to Polars DataFrame
    pairs_df = abutils.io.to_polars(
        pairs, 
    )

    # convert Polars DataFrame back to list of pairs
    pairs = abutils.io.from_polars(pairs_df)

|  

api
------------------


.. _to-pandas:  

.. autofunction:: abutils.io.to_pandas

.. _from-pandas:  

.. autofunction:: abutils.io.from_pandas

.. _to-polars:  

.. autofunction:: abutils.io.to_polars

.. _from-polars:  

.. autofunction:: abutils.io.from_polars





.. _AIRR-C: https://docs.airr-community.org/en/stable/datarep/rearrangements.html
.. _Polars: https://pola.rs/
.. _Pandas: https://pandas.pydata.org/
