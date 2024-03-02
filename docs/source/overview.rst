overview
========

Antibody repertoire sequencing is an increasingly important tool for detailed characterization 
of the immune response to infection and immunization. We built abutils to provide a cohesive 
set of tools designed for the specific challenges inherent in working with antibody repertoire data. 
The components in abutils were designed to be flexible: equally at home when used used interactively 
(in a Jupyter Notebook, for example) or when integrated into more complex programs and/or pipelines 
(such as abstar_, which is capable of annotating billions of antibody sequences).
  

core models
-----------

To represent antibody repertoire data at varying levels of granularity, abutils provides three core models: 
``Sequence``, ``Pair``, and ``Lineage``. These models are organized heirarchically -- a ``Lineage`` is composed of one 
or more ``Pair`` objects, and a ``Pair`` is composed of one or more ``Sequence`` objects -- with each model 
containing functionality appropriate for their respective level of granularity. 

  * ``Sequence``: model for representing a single antibody sequnce (either heavy or light chain).
    Provides a means to store and access abstar annotations. Includes common methods of sequence
    manipulation, including slicing, reverse-complement, translation, and conversion to FASTA format. 
    The ``Sequence`` object is used extensively throughout the ab[x] toolkit.  
  * ``Pair``: model for representing paired (heavy and light) antibody sequences. Comprised of one 
    or more ``Sequence`` objects. Heavily used in scab_, which is our toolkit for analyzing adaptive immune 
    single cell datasets.
  * ``Lineage``: model for representing an antibody clonal lineage. Comprised of one or more ``Pair``
    objects. Includes methods for lineage manipulation, including generating dot alignments and UCA calculation.


tools (abutils.tl)
------------------

In addition to the core models, abutils provides a number of commonly used functions. 
These functions are widely used throughput the ab[x] toolkit and are suitable for incorporation 
into custom pipelines or for use when performing interactive analyses:

  * :ref:`pairwise alignment <pairwise-alignment>`: local (Smith-Waterman), global (Needleman-Wunsch) and semi-global pairwise sequence alignment using parasail_
  
  * :ref:`multiple sequence alignment <multiple-sequence-alignment>` using MAFFT_ or MUSCLE_

  * :ref:`clustering <clustering>`: identity-based sequence clustering with VSEARCH_, CDHIT_, or MMseqs2_

  * :ref:`phylogeny <phylogeny>`: computing phylogenies with FastTree_ or IgPhyML_, tree drawing with baltic_



plots (abutils.pl)
------------------

``abutils`` provides a number of plotting functions for visualizing antibody repertoire data. These functions are
designed to be generalizable and flexible, and are suitable for use in both interactive and programmatic contexts.
Primary plotting functions include:

  * :ref:`bar <bar-plot>`: plot categorical data or frequency distributions as a bar plot

  * :ref:`scatter <scatter-plot>`: plot two-dimensional data as a scatter plot

  * :ref:`kde <pkde-plot>`: plot one- or two-dimensional data as a kernel density estimate

  * :ref:`donut <donut-plot>`: plot categorical data (such as lineages or germline genes) as a donut plot




.. _abstar: https://github.com/briney/abstar
.. _scab: https://github.com/briney/scab
.. _parasail: https://github.com/jeffdaily/parasail-python
.. _MAFFT: https://mafft.cbrc.jp/alignment/software/
.. _MUSCLE: https://www.drive5.com/muscle/
.. _VSEARCH: https://github.com/torognes/vsearch
.. _CDHIT: http://weizhongli-lab.org/cd-hit/
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _FastTree: http://www.microbesonline.org/fasttree/
.. _IgPhyML: https://github.com/kbhoehn/IgPhyML
.. _baltic: https://github.com/evogytis/baltic
