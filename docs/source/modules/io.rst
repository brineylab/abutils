.. _io:

sequence I/O
====================

``abutils`` provides a set of functions for :ref:`reading <read_sequences>` and :ref:`writing <write_sequences>` 
sequence data to and from various file formats. Additionally, we can :ref:`convert <convert_sequences>` lists of 
``Pair`` and ``Sequence`` objects to and from Pandas or Polars DataFrames.

|  

sequence annotations
-----------------------

``abutils`` follows the AIRR-C standard for sequence annotations. In tabular format, such as 
tab-delimited (the official AIRR format), CSV, or Parquet, sequence annotations appear as follows,
with one sequence per row:

.. csv-table::
   :header: "sequence_id", "sequence", "sequence_aa", "..."
   :widths: 20, 20, 20, 5
   :align: left

   "sequence1", "ATCG...", "EVQLVE...", "..."
   "sequence2", "ATCG...", "QVQLVE...", "..."
   "sequence3", "ATCG...", "EVQLVE...", "..."

|  

pair annotations
--------------------

``abutils`` uses a custom extentension of the AIRR-C standard for pair annotations. Each row contains 
a heavy/light chain pair. All AIRR-C fields are supported for each chain, with heavy chain annotation fields
appended with ``":0"`` and light chain annotation fields appended with ``":1"``. Additionally, a ``name`` field 
is included to allow for naming the pair independently of either sequence chain:

.. csv-table::
   :header: "name", "sequence_id:0", "sequence:0", "...", "sequence_id:1", "sequence:1", "..."
   :widths: 10, 20, 20, 5, 20, 20, 5
   :align: left

   "pair1", "sequence1_heavy", "ATGC...", "...", "sequence1_light", "ATGC...", "..."
   "pair2", "sequence2_heavy", "ATGC...", "...", "sequence2_light", "ATGC...", "..."
   "pair3", "sequence3_heavy", "ATGC...", "...", "sequence3_light", "ATGC...", "..."

|  

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
