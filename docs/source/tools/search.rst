.. _mmseqs-search:

search
=========

``abutils`` provides functions for searching sequences against a database of target sequences
using MMseqs2_. This allows for fast and efficient similarity search, which is useful for tasks
like sequence identification, annotation, and homology detection.

|  

.. csv-table:: 
   :header: "search method", "function"
   :widths: 10, 16

   "MMseqs2", :ref:`abutils.tools.mmseqs_search() <mmseqs-search>`

examples
---------

**search with MMseqs2**

The MMseqs2 search function can accept a path to a FASTA file, an MMseqs2 database, a :class:`abutils.Sequence` object,
or an iterable of :class:`abutils.Sequence` objects for both query and target sequences.

.. code-block:: python

    import abutils

    # search sequences against a target database
    results = abutils.tools.mmseqs_search(
        query='path/to/query_sequences.fasta',
        target='path/to/target_sequences.fasta',
        output_path='path/to/output.tsv'
    )

| 

**customize search parameters**

MMseqs2 search can be customized with various parameters to control sensitivity, format, and performance.

.. code-block:: python

    import abutils

    # search with customized parameters
    results = abutils.tools.mmseqs_search(
        query='path/to/query_sequences.fasta',
        target='path/to/target_sequences.fasta',
        output_path='path/to/output.tsv',
        search_type=1,  # amino acid search
        max_seqs=100,   # maximum hits per query
        max_evalue=1e-5,  # stricter E-value cutoff
        sensitivity=7.5,  # higher sensitivity
        format_mode=4,    # BLAST-TAB + column headers
        threads=8         # use 8 threads
    )

| 

**customizing output format**

You can customize the output format to include specific columns.

.. code-block:: python

    import abutils

    # customize the output format 
    results = abutils.tools.mmseqs_search(
        query='path/to/query_sequences.fasta',
        target='path/to/target_sequences.fasta',
        output_path='path/to/output.tsv',
        format_mode=4,  # BLAST-TAB + column headers
        format_output="query,target,evalue,pident,qcov,tcov" 
    )

| 

api
-------

.. _mmseqs-search:

.. autofunction:: abutils.tools.mmseqs_search

.. _MMseqs2:
    https://github.com/soedinglab/MMseqs2 