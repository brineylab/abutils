.. _preprocessing:

preprocessing
=============

``abutils`` provides several functions for preprocessing sequence data, including merging paired-end
FASTQ files. The preprocessing module can handle different file naming schemas (Illumina and Element)
and supports various merging algorithms. 

The primary function is ``abutils.tools.merge_fastqs()``, which handles the process of organizing, 
grouping, and merging paired-end FASTQ files.

|  

.. csv-table:: 
   :header: "preprocessing method", "function"
   :widths: 10, 16

   "Merge paired-end reads", :ref:`abutils.tools.merge_fastqs() <merge-fastqs>`
   "Merge with fastp", :ref:`abutils.tools.merge_fastqs_fastp() <merge-fastqs-fastp>`
   "Merge with vsearch", :ref:`abutils.tools.merge_fastqs_vsearch() <merge-fastqs-vsearch>`

examples
---------

**merge paired-end FASTQ files using fastp**

By default, ``merge_fastqs()`` uses fastp to merge paired-end reads, which also performs
quality filtering and adapter trimming.

.. code-block:: python

    import abutils

    # merge paired-end FASTQ files from a directory
    merged_files = abutils.tools.merge_fastqs(
        files='path/to/fastq/directory',
        output_directory='path/to/output',
        schema='illumina',  # file naming schema: 'illumina' or 'element'
        compress_output=True,  # output gzipped FASTQ files
        verbose=True
    )

| 

**customize quality trimming and adapter removal**

You can customize the quality trimming and adapter removal parameters.

.. code-block:: python

    import abutils

    # merge with customized quality trimming and adapter removal
    merged_files = abutils.tools.merge_fastqs(
        files='path/to/fastq/directory',
        output_directory='path/to/output',
        trim_adapters=True,
        adapter_file='path/to/adapters.fasta',  # custom adapter sequences
        quality_trim=True,
        window_size=5,  # sliding window size
        quality_cutoff=15,  # quality threshold
        minimum_overlap=20,  # minimum overlap between reads
        log_directory='path/to/logs'  # save fastp reports
    )

| 

**direct low-level merging with specific algorithms**

For more control, you can directly use the algorithm-specific merge functions. This requires
specifying each of the input files and output file paths, rather than simply input and output directories.

.. code-block:: python

    import abutils

    # merge a specific pair of files with fastp
    abutils.tools.merge_fastqs_fastp(
        forward='path/to/sample_R1.fastq.gz',
        reverse='path/to/sample_R2.fastq.gz',
        merged='path/to/output/sample.fastq.gz',
        minimum_overlap=25,
        allowed_mismatches=3,
        trim_adapters=True,
        quality_trim=True
    )

    # or with vsearch
    abutils.tools.merge_fastqs_vsearch(
        forward='path/to/sample_R1.fastq.gz',
        reverse='path/to/sample_R2.fastq.gz',
        merged_file='path/to/output/sample.fastq',
        output_format='fastq',
        minimum_overlap=25,
        allowed_mismatches=3
    )

| 

api
-------

.. _merge-fastqs:

.. autofunction:: abutils.tools.merge_fastqs

.. _merge-fastqs-fastp:

.. autofunction:: abutils.tools.merge_fastqs_fastp

.. _merge-fastqs-vsearch:

.. autofunction:: abutils.tools.merge_fastqs_vsearch 