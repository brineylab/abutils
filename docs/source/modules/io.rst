.. _io:

sequence I/O
====================

``abutils`` provides a set of functions for reading and writing sequence data to and from various file formats. Additionally,
we can convert lists of ``Pair`` and ``Sequence`` objects to and from Pandas or Polars DataFrames.


read
----------------------------

.. table::
   :class: centered-table

   +----------+--------------------------------------+----------------------------------------+
   | format   | function                             | notes                                  |
   +==========+======================================+========================================+
   | FASTA/Q  | :ref:`read_fastx() <read-fastx>`     | supports gzip-compressed files         |
   |          |                                      |                                        |
   +          +--------------------------------------+----------------------------------------+
   |          | :ref:`parse_fastx() <parse-fastx>`   | supports gzip-compressed files         |
   |          |                                      |                                        |
   +----------+--------------------------------------+----------------------------------------+
   | FASTA    | :ref:`read_fasta() <read-fasta>`     | supports gzip-compressed files         |
   |          |                                      |                                        |
   +          +--------------------------------------+----------------------------------------+
   |          | :ref:`parse_fasta() <parse-fasta>`   | supports gzip-compressed files         |
   |          |                                      |                                        |
   +----------+--------------------------------------+----------------------------------------+
   | FASTQ    | :ref:`read_fastq() <read-fastq>`     | supports gzip-compressed files         |
   |          |                                      |                                        |
   +          +--------------------------------------+----------------------------------------+
   |          | :ref:`parse_fastq() <parse-fastq>`   | supports gzip-compressed files         |
   |          |                                      |                                        |
   +----------+--------------------------------------+----------------------------------------+
   | AIRR     | :ref:`read_airr() <read-airr>`       | supports ``Sequence`` objects only     |
   |          |                                      |                                        |
   +----------+--------------------------------------+----------------------------------------+
   | Parquet  | :ref:`read_parquet() <read-parquet>` | supports ``Sequence`` or ``Pair``      |
   |          |                                      | objects                                |
   +----------+--------------------------------------+----------------------------------------+
   | CSV      | :ref:`read_csv() <read-csv>`         | supports ``Sequence`` or ``Pair``      |
   |          |                                      | objects                                |
   +----------+--------------------------------------+----------------------------------------+


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

.. table::
   :class: centered-table

   +---------+--------------------------------------+----------------------------------------+
   | format  | function                             | notes                                  |
   +=========+======================================+========================================+
   | Pandas  | :ref:`to_pandas()                    | supports ``Sequence`` or ``Pair``      |
   |         | <convert_sequences#to_pandas>`       | objects                                |
   +         +--------------------------------------+----------------------------------------+
   |         | :ref:`from_pandas()                  | supports ``Sequence`` or ``Pair``      |
   |         | <convert_sequences#from_pandas>`     | objects                                |
   +---------+--------------------------------------+----------------------------------------+
   | Polars  | :ref:`to_polars()                    | supports ``Sequence`` or ``Pair``      |
   |         | <convert_sequences#to_polars>`       | objects                                |
   +         +--------------------------------------+----------------------------------------+
   |         | :ref:`from_polars()                  | supports ``Sequence`` or ``Pair``      |
   |         | <convert_sequences#from_polars>`     | objects                                |
   +---------+--------------------------------------+----------------------------------------+



.. toctree::
    :hidden:

    read <read_sequences>
    write <write_sequences>
    convert <convert_sequences>
