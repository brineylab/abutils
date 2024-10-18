

read/parse sequences
==============================


.. csv-table:: 
   :header: "format", "function", "notes"
   :widths: 15, 10, 30

   "FASTA/Q", :ref:`read_fastx() <read-fastx>`, "also supports gzipped files"
   "FASTA", :ref:`read_fasta() <read-fasta>`, "also supports gzipped files"
   "FASTQ", :ref:`read_fastq() <read-fastq>`, "also supports gzipped files"
   "AIRR", :ref:`read_airr() <read-airr>`, "only supports ``Sequence`` objects"
   "Parquet", :ref:`read_parquet() <read-parquet>`, "supports ``Sequence`` or ``Pair`` objects"
   "CSV", :ref:`read_csv() <read-csv>`, "supports ``Sequence`` or ``Pair`` objects"



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


