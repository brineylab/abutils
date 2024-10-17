

read sequences
==============================


.. csv-table:: read
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA/Q", "`read_fastx() <read_fastx>`_", "also supports gzipped files"
   "FASTA", "`read_fasta() <read_fasta>`_", "also supports gzipped files"
   "FASTQ", "`read_fastq() <read_fastq>`_", "also supports gzipped files"
   "Parquet", "`read_parquet() <read_parquet>`_", "supports ``Sequence`` or ``Pair`` objects"
   "CSV", "`read_csv() <read_csv>`_", "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", "`read_airr() <read_airr>`_", "only supports ``Sequence`` objects"



.. _read_fastx:
.. autofunction:: abutils.core.sequence.read_fastx

.. _read_fasta:
.. autofunction:: abutils.core.sequence.read_fasta

.. _read_fastq:
.. autofunction:: abutils.core.sequence.read_fastq

.. _read_airr:
.. autofunction:: abutils.core.sequence.read_airr

.. _read_csv:
.. autofunction:: abutils.core.sequence.read_csv

.. _read_parquet:
.. autofunction:: abutils.core.sequence.read_parquet



