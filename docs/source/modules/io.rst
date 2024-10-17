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

   "FASTA/Q", :ref:`read_fastx() <read-fastx>`, "also supports gzipped files"
   "FASTA", :ref:`read_fasta() <read-fasta>`, "also supports gzipped files"
   "FASTQ", :ref:`read_fastq() <read-fastq>`, "also supports gzipped files"
   "Parquet", :ref:`read_parquet() <read-parquet>`, "supports ``Sequence`` or ``Pair`` objects"
   "CSV", :ref:`read_csv() <read-csv>`, "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", :ref:`read_airr() <read-airr>`, "only supports ``Sequence`` objects"


write
----------------------------

.. csv-table::
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA", :ref:`to_fasta() <to-fasta>`, "supports ``Sequence`` or ``Pair`` objects"
   "FASTQ", :ref:`to_fastq() <to-fastq>`, "supports ``Sequence`` or ``Pair`` objects"
   "Parquet", :ref:`to_parquet() <to-parquet>`, "supports ``Sequence`` or ``Pair`` objects"
   "CSV", :ref:`to_csv() <to-csv>`, "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", :ref:`to_airr() <to-airr>`, "only supports ``Sequence`` objects"


convert
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 30

   * - format
     - function
     - notes
   * - Pandas
     - :ref:`to_pandas() <convert_sequences#to_pandas>`
     - supports ``Sequence`` or ``Pair`` objects
   * - 
     - :ref:`from_pandas() <convert_sequences#from_pandas>`
     - supports ``Sequence`` or ``Pair`` objects
   * - Polars
     - :ref:`to_polars() <convert_sequences#to_polars>`
     - supports ``Sequence`` or ``Pair`` objects
   * - 
     - :ref:`from_polars() <convert_sequences#from_polars>`
     - supports ``Sequence`` or ``Pair`` objects


.. toctree::
    :hidden:

    read <read_sequences>
    write <write_sequences>
    convert <convert_sequences>
