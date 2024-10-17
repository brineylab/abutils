

read sequences
==============================


.. csv-table:: read
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA/Q", :ref:`read_fastx() <read_fastx>`, "also supports gzipped files"
   "FASTA", :ref:`read_fasta() <read_fasta>`, "also supports gzipped files"
   "FASTQ", :ref:`read_fastq() <read_fastq>`, "also supports gzipped files"
   "Parquet", :ref:`read_parquet() <read_parquet>`, "supports ``Sequence`` or ``Pair`` objects"
   "CSV", :ref:`read_csv() <read_csv>`, "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", :ref:`read_airr() <read_airr>`, "only supports ``Sequence`` objects"



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



