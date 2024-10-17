.. _io:

sequence I/O
====================

``abutils`` provides a set of functions for reading and writing sequence data to and from various file formats. Additionally,
we can convert lists of ``Pair`` and ``Sequence`` objects to and from Pandas or Polars DataFrames.


read
----------------------------

.. csv-table::
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA/Q", "`read_fastx() <read-fastx>`_", "also supports gzipped files"
   "FASTA", "`read_fasta() <read-fasta>`_", "also supports gzipped files"
   "FASTQ", "`read_fastq() <read-fastq>`_", "also supports gzipped files"
   "Parquet", "`read_parquet() <read-parquet>`_", "supports ``Sequence`` or ``Pair`` objects"
   "CSV", "`read_csv() <read-csv>`_", "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", "`read_airr() <read-airr>`_", "only supports ``Sequence`` objects"


write
----------------------------

.. csv-table::
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA", "`to_fasta() <to-fasta>`_", "supports ``Sequence`` or ``Pair`` objects"
   "FASTQ", "`to_fastq() <to-fastq>`_", "supports ``Sequence`` or ``Pair`` objects"
   "Parquet", "`to_parquet() <to-parquet>`_", "supports ``Sequence`` or ``Pair`` objects"
   "CSV", "`to_csv() <to-csv>`_", "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", "`to_airr() <to-airr>`_", "only supports ``Sequence`` objects"


convert
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 30

   * - format
     - function
     - notes
   * - Pandas
     - `to_pandas() <convert_sequences#to_pandas>`_
     - supports ``Sequence`` or ``Pair`` objects
   * - 
     - `from_pandas() <convert_sequences#from_pandas>`_
     - supports ``Sequence`` or ``Pair`` objects
   * - Polars
     - `to_polars() <convert_sequences#to_polars>`_
     - supports ``Sequence`` or ``Pair`` objects
   * - 
     - `from_polars() <convert_sequences#from_polars>`_
     - supports ``Sequence`` or ``Pair`` objects


.. toctree::
    :hidden:

    read <read_sequences>
    write <write_sequences>
    convert <convert_sequences>
