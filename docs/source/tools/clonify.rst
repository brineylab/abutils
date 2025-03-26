.. _clonify:

clonify
=========

``abutils`` provides functions for assigning antibody sequences to clonal lineages using the clonify_ [Briney16]_ 
algorithm. This algorithm uses a combination of CDR3 sequence similarity and shared somatic hypermutation 
patterns to group sequences into B cell clonal lineages.

The primary function is ``abutils.tl.clonify()``, which handles lineage assignment at scale, with support
for different input/output formats and parallel processing.

|  

.. csv-table:: 
   :header: "lineage assignment method", "function"
   :widths: 10, 16

   "Clonify algorithm", :ref:`abutils.tl.clonify() <clonify-function>`
   "Pairwise distance calculation", :ref:`abutils.tl.pairwise_distance() <pairwise-distance-function>`

examples
---------

**basic lineage assignment**

``clonify()`` can accept a variety of input formats, including paths to AIRR-formatted TSV files, 
Parquet files, or lists of :class:`abutils.Sequence` objects.

.. code-block:: python

    import abutils

    # clonal assignment using default parameters
    lineages = abutils.tl.clonify(
        sequences='path/to/airr_data.tsv',
        output_path='path/to/output_with_lineages.tsv',
        verbose=True
    )

| 

**customizing lineage assignment parameters**

You can customize the parameters that control lineage assignment sensitivity and specificity.

.. code-block:: python

    import abutils

    # customize lineage assignment parameters
    lineages = abutils.tl.clonify(
        sequences='path/to/airr_data.tsv',
        output_path='path/to/output_with_lineages.tsv',
        distance_cutoff=0.32,             # stricter distance threshold
        shared_mutation_bonus=0.4,        # increased bonus for shared mutations
        length_penalty_multiplier=2.5,    # increased penalty for CDR3 length differences
        group_by_v=True,                  # group by V-gene before assignment
        group_by_j=True,                  # group by J-gene before assignment
        verbose=True
    )

| 

**working with paired heavy and light chain data**

For paired data, you can use light chain information in the lineage assignment process.

.. code-block:: python

    import abutils

    # lineage assignment with paired heavy/light chain data
    lineages = abutils.tl.clonify(
        sequences='path/to/paired_data.parquet',
        output_path='path/to/output_with_lineages.parquet',
        output_fmt='parquet',
        group_by_light_chain_vj=True,     # also group by light chain V/J genes
        n_processes=8                     # use 8 processes for parallel computation
    )

| 

api
-------

.. _clonify-function:

.. autofunction:: abutils.tools.clonify

.. _pairwise-distance-function:

.. autofunction:: abutils.tools.pairwise_distance

.. _clonify:
    https://doi.org/10.1038/srep23901

.. [Briney16] Bryan Briney, Khoa Le, Jiang Zhu, and Dennis R Burton. Clonify: unseeded antibody lineage assignment from next-generation sequencing data. *Scientific Reports* 2016. https://doi.org/10.1038/srep23901 