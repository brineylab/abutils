.. abutils documentation master file, created by
   sphinx-quickstart on Fri Jul 27 10:54:47 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

abutils: utilities for AIRR analysis
====================================================================


Antibody repertoire sequencing is an increasingly important tool for detailed characterization 
of the immune response to infection and immunization. We built abutils to provide a cohesive 
set of tools designed for the specific challenges inherent in working with antibody repertoire data. 
The components in abutils were designed to be flexible: equally at home when used used interactively 
(in a Jupyter Notebook, for example) or when integrated into more complex programs and/or pipelines 
(such as abstar_, which is capable of annotating billions of antibody sequences).

| 

core models
----------------------

To represent antibody repertoire data at varying levels of granularity, ``abutils`` provides three core models:  

  * :ref:`Sequence <sequence-class>`: model for representing a single antibody sequnce (either heavy or light chain).
    Provides a means to store and access abstar annotations. Includes common methods of sequence
    manipulation, including slicing, reverse-complement, and conversion to FASTA format. The ``Sequence``
    object is used extensively throughout the ab[x] toolkit.
  * :ref:`Pair <pair-class>`: model for representing paired (heavy and light) antibody sequences. Comprised of one 
    or more ``Sequence`` objects. Heavily used in scab_, which is our toolkit for analyzing adaptive immune 
    single cell datasets.
  * :ref:`Lineage <lineage-class>`: model for representing an antibody clonal lineage. Comprised of one or more ``Pair``
    objects. Includes methods for lineage manipulation, including generating dot alignments and UCA calculation.

 
These models are heirarchical -- a ``Lineage`` is composed of one 
or more ``Pair`` objects, a ``Pair`` is composed of one or more ``Sequence`` objects -- and contain methods 
appropriate for each level of granularity.

| 

data (``abutils.io``)
--------------------------------

To simplify data manipulation and facilitate the integration of ``abutils`` into 
existing pipelines, ``abutils`` provides a set of functions for reading and writing 
sequence data:

  * :ref:`read <read-sequences>`: read sequences from a variety of file formats, including FASTA, FASTQ, AIRR-C, and others.
  
  * :ref:`write <write-sequences>`: write sequences to a variety of file formats, including FASTA, FASTQ, AIRR-C, and others.
  
  * :ref:`convert <convert-sequences>`: convert between ``Sequence`` or ``Pair`` objects and Pandas_ or Polars_ DataFrames.

  * :ref:`paths <path>`: functions for working with file paths and directories.

All of the IO functions are accessible via ``abutils.io``.

| 

tools (``abutils.tl``)
--------------------

In addition to the core models, abutils provides a number of commonly used functions. 
These functions are widely used throughput the ab[x] toolkit and can be easily integrated
into custom pipelines or for use when performing interactive analyses:

  * :ref:`pairwise alignment <pairwise-alignment>`: local (Smith-Waterman), global (Needleman-Wunsch) and semi-global pairwise sequence alignment using parasail_. 
  
  * :ref:`multiple sequence alignment <multiple-sequence-alignment>` using MAFFT_, MUSCLE_, or FAMSA_

  * :ref:`clustering <clustering>`: identity-based sequence clustering with VSEARCH_, CDHIT_, or MMseqs2_

  * :ref:`search <search>`: fast sequence similarity search with MMseqs2_
  
  * :ref:`preprocessing <preprocessing>`: preprocessing functions for sequence data, including merging paired-end reads
  
  * :ref:`clonify <clonify>`: assigning antibody sequences to clonal lineages using the clonify_ algorithm

  * :ref:`phylogeny <phylogeny>`: computing phylogenies with FastTree_ or IgPhyML_, tree drawing with baltic_

All of the tool functions are accessible via ``abutils.tl``.

| 

plots (``abutils.pl``)
--------------------

``abutils`` provides a number of plotting functions for visualizing antibody repertoire data.
These functions are built on top of matplotlib and seaborn and are designed to be easily
integrated into custom analyses or pipelines. Plotting functions are designed to work with
``Sequence``, ``Pair``, and ``Lineage`` objects, and fully support AIRR-C annotation formats for
plotting adaptive immune receptor features like CDR3 length distributions and germline gene usage.

  * :ref:`bar <bar-plot>`: plot categorical data or frequency distributions as a bar plot

  * :ref:`scatter <scatter-plot>`: plot two-dimensional data as a scatter plot

  * :ref:`kde <kde-plot>`: plot one- or two-dimensional data as a kernel density estimate

  * :ref:`donut <donut-plot>`: plot categorical data (such as lineages or germline genes) as a donut plot

All of the plotting functions are accessible via ``abutils.pl``



colors (``abutils.cl``)
-------------------------

``abutils`` provides utility functions that are generally useful when working with antibody repertoire data.
These include functions for monitoring multiprocessing jobs, creating and manipulating color palettes, and more.

  * :ref:`colors <colors>`: comprehensive functions for working with colors and color palettes, accessible through the ``abutils.cl`` module



.. _abstar: https://github.com/briney/abstar
.. _scab: https://github.com/briney/scab
.. _parasail: https://github.com/jeffdaily/parasail-python
.. _MAFFT: https://mafft.cbrc.jp/alignment/software/
.. _MUSCLE: https://www.drive5.com/muscle/
.. _FAMSA: https://github.com/MikkelSchubert/FAMSA
.. _VSEARCH: https://github.com/torognes/vsearch
.. _CDHIT: http://weizhongli-lab.org/cd-hit/
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _FastTree: http://www.microbesonline.org/fasttree/
.. _IgPhyML: https://github.com/kbhoehn/IgPhyML
.. _baltic: https://github.com/evogytis/baltic
.. _Pandas: https://pandas.pydata.org/
.. _Polars: https://pola.rs/
.. _clonify: https://doi.org/10.1038/srep23901



|

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: getting started
   
   installation


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: models

   modules/sequence
   modules/pair


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: data

   io <modules/io>
   path <modules/path>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: tools

   tools/alignment
   tools/clonify
   tools/clustering
   tools/phylogeny
   tools/preprocessing
   tools/search

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: plots

   plots/bar
   plots/scatter
   plots/kde
   plots/donut

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: colors
  
    modules/color

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: about

   license


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: related projects

   abstar <https://github.com/brineylab/abstar>
   scab <https://github.com/brineylab/scab>


index
-----

* :ref:`genindex`
* :ref:`modindex`
