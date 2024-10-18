

write sequence data
=========================

``abutils`` provides functions for writing sequence data to a variety of commonly 
used file formats. This includes raw sequence data in FASTA or FASTQ format as well as 
annotated sequence data in AIRR-C_, CSV, or Parquet formats.

|  


.. csv-table:: 
   :header: "format", "function", "notes"
   :align: left
   :widths: 10, 12, 24
   :width: 100%

   "FASTA", :ref:`to_fasta() <to-fasta>`, "supports ``Sequence`` or ``Pair`` objects"
   "FASTQ", :ref:`to_fastq() <to-fastq>`, "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", :ref:`to_airr() <to-airr>`, "only supports ``Sequence`` objects"
   "Parquet", :ref:`to_parquet() <to-parquet>`, "supports ``Sequence`` or ``Pair`` objects"
   "CSV", :ref:`to_csv() <to-csv>`, "supports ``Sequence`` or ``Pair`` objects"

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

.. _to-fasta:  

.. autofunction:: abutils.io.to_fasta

.. _to-fastq:  

.. autofunction:: abutils.io.to_fastq

.. _to-airr:  

.. autofunction:: abutils.io.to_airr

.. _to-parquet:  

.. autofunction:: abutils.io.to_parquet

.. _to-csv:

.. autofunction:: abutils.io.to_csv

