

write sequences
=========================




.. csv-table:: 
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA", :ref:`to_fasta() <to-fasta>`, "supports ``Sequence`` or ``Pair`` objects"
   "FASTQ", :ref:`to_fastq() <to-fastq>`, "supports ``Sequence`` or ``Pair`` objects"
   "AIRR", :ref:`to_airr() <to-airr>`, "only supports ``Sequence`` objects"
   "Parquet", :ref:`to_parquet() <to-parquet>`, "supports ``Sequence`` or ``Pair`` objects"
   "CSV", :ref:`to_csv() <to-csv>`, "supports ``Sequence`` or ``Pair`` objects"


.. _to-fasta:  

.. autofunction:: abutils.core.sequence.to_fasta

.. _to-fastq:  

.. autofunction:: abutils.core.sequence.to_fastq

.. _to-airr:  

.. autofunction:: abutils.io.to_airr

.. _to-parquet:  

.. autofunction:: abutils.io.to_parquet

.. _to-csv:

.. autofunction:: abutils.io.to_csv

