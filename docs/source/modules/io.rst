.. _io:

sequence I/O
====================

``abutils`` provides a set of functions for reading and writing sequence data to and from various file formats. Additionally,
we can convert lists of ``Pair`` and ``Sequence`` objects to and from Pandas or Polars DataFrames.


read
----------------------------

.. table::
   :align: left
   :widths: 10, 12, 24
   :width: 100%

   +----------+--------------------------------------+----------------------------------------+
   | format   | function                             | notes                                  |
   +==========+======================================+========================================+
   | FASTA/Q  | :ref:`read_fastx() <read-fastx>`     | returns a list of ``Sequence`` objects |
   |          |                                      |                                        |
   +          +--------------------------------------+----------------------------------------+
   |          | :ref:`parse_fastx() <parse-fastx>`   | yields single ``Sequence`` objects     |
   |          |                                      |                                        |
   +----------+--------------------------------------+----------------------------------------+
   | FASTA    | :ref:`read_fasta() <read-fasta>`     | returns a list of ``Sequence`` objects |
   |          |                                      |                                        |
   +          +--------------------------------------+----------------------------------------+
   |          | :ref:`parse_fasta() <parse-fasta>`   | yields single ``Sequence`` objects     |
   |          |                                      |                                        |
   +----------+--------------------------------------+----------------------------------------+
   | FASTQ    | :ref:`read_fastq() <read-fastq>`     | returns a list of ``Sequence`` objects |
   |          |                                      |                                        |
   +          +--------------------------------------+----------------------------------------+
   |          | :ref:`parse_fastq() <parse-fastq>`   | yields single ``Sequence`` objects     |
   |          |                                      |                                        |
   +----------+--------------------------------------+----------------------------------------+
   | AIRR     | :ref:`read_airr() <read-airr>`       | only supports ``Sequence`` objects     |
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
   :widths: 10, 12, 24
   :width: 100%
   :align: left

   "FASTA", :ref:`to_fasta() <to-fasta>`, "supports ``Sequence`` or ``Pair`` objects"
   "FASTQ", :ref:`to_fastq() <to-fastq>`, "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", :ref:`to_airr() <to-airr>`, "only supports ``Sequence`` objects"
   "Parquet", :ref:`to_parquet() <to-parquet>`, "supports ``Sequence`` or ``Pair`` objects"
   "CSV", :ref:`to_csv() <to-csv>`, "supports ``Sequence`` or ``Pair`` objects"


convert
----------------------------

.. table::
   :align: left
   :widths: 10, 12, 24
   :width: 100%

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
