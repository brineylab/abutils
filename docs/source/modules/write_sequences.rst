

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

fasta/q files
------------------

``abutils`` can write lists of ``Seuqence`` or ``Pair`` objects to FASTA or FASTQ files:

.. warning::

    While ``abutils`` can write both ``Sequence`` and ``Pair`` objects to various formats, the input must contains
    only one type of object. For example, you cannot mix ``Sequence`` and ``Pair`` objects in the 
    same list.

.. code-block:: python

    # write list of sequences to FASTA file
    abutils.io.to_fasta(
        sequences, 
        "my-output-file.fasta"
    )

    # write list of pairs to FASTQ file
    abutils.io.to_fastq(
        pairs, 
        "my-paired-output-file.fastq"
    )

|  

annotated sequence files
------------------------------
``to_airr()`` can write ``Sequence`` objects to AIRR-C_ formatted (tab-delimited) files:

.. code-block:: python

    # write list of sequences to AIRR-C file
    abutils.io.to_airr(
        sequences, 
        "my-airr-output-file.tsv"
    )

|  

``to_parquet()`` can write ``Sequence`` or ``Pair`` objects to Parquet files:

.. code-block:: python

    # write list of sequences to Parquet file
    abutils.io.to_parquet(
        sequences, 
        "my-parquet-output-file.parquet"
    )

    # write list of pairs to Parquet file
    abutils.io.to_parquet(
        pairs, 
        "my-paired-parquet-output-file.parquet"
    )


|  

``to_csv()`` can write ``Sequence`` or ``Pair`` objects to CSV files:

.. code-block:: python

    # write list of sequences to CSV file
    abutils.io.to_csv(
        sequences, 
        "my-csv-output-file.csv"

    # write list of pairs to CSV file
    abutils.io.to_csv(
        pairs, 
        "my-paired-csv-output-file.csv"
    )


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

