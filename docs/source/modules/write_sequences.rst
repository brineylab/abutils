

write sequences
=========================


.. csv-table:: write
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA", "`to_fasta() <to-fasta>`_", "supports ``Sequence`` or ``Pair`` objects"
   "FASTQ", "`to_fastq() <to-fastq>`_", "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", "`to_airr() <to-airr>`_", "only supports ``Sequence`` objects"
   "Parquet", "`to_parquet() <to-parquet>`_", "supports ``Sequence`` or ``Pair`` objects"
   "CSV", "`to_csv() <to-csv>`_", "supports ``Sequence`` or ``Pair`` objects"


.. _to-fasta:
.. autofunction:: abutils.core.sequence.to_fasta

.. _to-fastq:
.. autofunction:: abutils.core.sequence.to_fastq

.. _to-airr:
.. autofunction:: abutils.core.sequence.to_airr

.. _to-parquet:
.. autofunction:: abutils.core.sequence.to_parquet

.. _to-csv:
.. autofunction:: abutils.core.sequence.to_csv

