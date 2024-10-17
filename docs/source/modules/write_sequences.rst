

write sequences
=========================


.. csv-table:: write
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA", "`to_fasta() <to_fasta>`_", "supports ``Sequence`` or ``Pair`` objects"
   "FASTQ", "`to_fastq() <to_fastq>`_", "supports ``Sequence`` or ``Pair`` objects"
   "Parquet", "`to_parquet() <to_parquet>`_", "supports ``Sequence`` or ``Pair`` objects"
   "CSV", "`to_csv() <to_csv>`_", "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", "`to_airr() <to_airr>`_", "only supports ``Sequence`` objects"


.. _to_fasta:
.. autofunction:: abutils.core.sequence.to_fasta

.. _to_fastq:
.. autofunction:: abutils.core.sequence.to_fastq

.. _to_parquet:
.. autofunction:: abutils.core.sequence.to_parquet

.. _to_csv:
.. autofunction:: abutils.core.sequence.to_csv

.. _to_airr:
.. autofunction:: abutils.core.sequence.to_airr

