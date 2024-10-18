

read/parse sequences
==============================

``abutils`` provides functions for reading/parsing sequence data from a variety of commonly 
used file formats. The primary differences between ``read`` and ``parse`` functions are:

- ``read`` functions read an entire file into memory and return a list of ``Sequence`` objects.
- ``parse`` functions yield ``Sequence`` objects one at a time.

Thus, ``parse`` functions are generally more memory efficient for large files, but ``read`` 
functions may be more convenient for smaller files or quick prototyping.

.. note::

    ``read_fastx()`` and ``parse_fastx()`` are the most flexible and can read/parse either
    FASTA or FASTQ files. This is particularly useful when building pipelines in which users 
    may want to process both file types or when the source file may not be known in advance.

.. code-block:: python

    # read entire file into memory
    sequences = abutils.io.read_fastx("sequences.fastq")

    # parse file one record at a time
    for sequence in abutils.io.parse_fastq("sequences.fastq"):
        print(sequence)




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



.. _read-fastx:  

.. autofunction:: abutils.io.read_fastx

.. _read-fasta:  

.. autofunction:: abutils.io.read_fasta

.. _read-fastq:  

.. autofunction:: abutils.io.read_fastq

.. _read-airr:  

.. autofunction:: abutils.io.read_airr

.. _read-parquet:  

.. autofunction:: abutils.io.read_parquet

.. _read-csv:  

.. autofunction:: abutils.io.read_csv


