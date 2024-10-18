

read sequence data
==============================


``abutils`` provides functions for reading/parsing sequence data from a variety of commonly 
used file formats. This includes raw sequence data in FASTA or FASTQ format as well as 
annotated sequence data in AIRR-C_, CSV, or Parquet formats.

|  

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


fasta/q files
------------------

The primary differences between ``read`` and ``parse`` functions are:

- ``read`` functions read an entire file into memory and return a list of ``Sequence`` objects.
- ``parse`` functions yield ``Sequence`` objects one at a time.

``parse`` functions are generally more memory efficient for large files, but ``read`` 
functions may be more convenient for smaller files or when sequences need to be processed
as a group rather than one at a time:

.. code-block:: python

    # read entire file into memory
    sequences = abutils.io.read_fasta("sequences.fasta")

    # parse file one record at a time
    for sequence in abutils.io.parse_fastq("sequences.fastq"):
        print(sequence)

|  

``read_fastx()`` and ``parse_fastx()`` are the most flexible and can read/parse either
FASTA or FASTQ files. This is particularly useful when building pipelines in which users 
may want to process both file types or when the source file may not be known in advance:

.. code-block:: python

    # FASTA file
    sequences = abutils.io.read_fastx("sequences.fasta")

    # FASTQ file
    for sequence in abutils.io.parse_fastx("sequences.fastq"):
        print(sequence)

|  

All of the FASTA/Q/X ``read`` and ``parse`` functions can handle gzip-compressed files automatically:

.. code-block:: python

    # FASTA file
    sequences = abutils.io.read_fastx("sequences.fasta.gz")

    # FASTQ file
    for sequence in abutils.io.parse_fastx("sequences.fastq.gz"):
        print(sequence)

|  

annotated sequence files
---------------------------

``read_airr()`` can read AIRR-C formatted sequence data from a tab-delimited file, 
returing a list of ``Sequence`` objects:

.. code-block:: python

    sequences = abutils.io.read_airr("sequences.tsv")

|  

``read_parquet()`` and ``read_csv()`` can read Parquet and CSV formatted annotated sequence data,
and expect the annotations to be in AIRR-C format -- the only difference is in the file format,
which can be either Parquet or CSV instead of the AIRR-C tab-delimited format: 

.. code-block:: python

    # read CSV file of annotated sequences
    sequences = abutils.io.read_csv("sequences.csv")

    # read Parquet file of annotated paired sequences
    pairs = abutils.io.read_parquet("pairs.parquet")

|  

Both ``read_csv()`` and ``read_parquet()`` support reading annotations from paired sequences, 
which is a custom extension of the AIRR-C format. Each row in the CSV or Parquet file 
contains annotations for both heavy and light chains. All annotation 
fields in the AIRR-C format are conserved for each chain, with heavy chains appending ``":0"`` 
to the end of each annotation field name and light chains appending ``":1"``. The row also contains
a ``"name"`` field so that the name of he paired sequence can be distinct from the names of the 
individual chains.

.. note::

    ``read_parquet()`` and ``read_csv()`` will automatically detect whether the input file
    contains ``Sequence`` or ``Pair`` objects based on the file schema.


.. code-block:: python

    # read CSV file of annotated paired sequences
    pairs = abutils.io.read_csv("pairs.csv")

    # read Parquet file of annotated paired sequences
    pairs = abutils.io.read_parquet("pairs.parquet")

|  

All of the functions for reading annotated sequence data include a ``match`` parameter that 
can be used to filter the sequences or pairs that are read from the file. This is useful 
when only a fraction of the sequences or pairs in the file are desired:

.. code-block:: python

    # read an AIRR file of sequences and return only those that use IGHV1-2
    sequences = abutils.io.read_airr(
        "sequences.tsv", 
        match={"v_gene": "IGHV1-2"},
    )

    # read Parquet file of paired sequences and return only those 
    # that have a productive heavy chain and light chain
    pairs = abutils.io.read_parquet(
        "pairs.parquet",
        match={
            "productive:0": True, 
            "productive:1": True,
        },
    )

|  

api
------------------


.. _read-fastx:  

.. autofunction:: abutils.io.read_fastx

.. _parse-fastx:  

.. autofunction:: abutils.io.parse_fastx

.. _read-fasta:  

.. autofunction:: abutils.io.read_fasta

.. _parse-fasta:  

.. autofunction:: abutils.io.parse_fasta

.. _read-fastq:  

.. autofunction:: abutils.io.read_fastq

.. _parse-fastq:  

.. autofunction:: abutils.io.parse_fastq

.. _read-airr:  

.. autofunction:: abutils.io.read_airr

.. _read-parquet:  

.. autofunction:: abutils.io.read_parquet

.. _read-csv:  

.. autofunction:: abutils.io.read_csv





.. _AIRR-C: https://docs.airr-community.org/en/stable/datarep/rearrangements.html