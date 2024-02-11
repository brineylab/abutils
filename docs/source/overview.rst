overview
========

With technical breakthroughs in the throughput and read-length of 
next-generation sequencing platforms, antibody repertoire sequencing 
is becoming an increasingly important tool for detailed characterization 
of the immune response to infection and immunization. Accordingly, 
there is a need for open, scalable software for the genetic analysis of 
repertoire-scale antibody sequence data.

We built ``abutils`` to provide a cohesive set of tools designed for the specific challenges inherent in 
working with antibody repertoire data. The components in ``abutils`` were designed to be flexible: 
equally at home when used used interactively (in a Jupyter Notebook, for example) or when 
integrated into more complex programs and/or pipelines (such as abstar_, which is capable of annotating 
billions of antibody sequences).

core models
-----------

To represent antibody repertoire data at varying levels of granularity, abutils provides three core models: 
``Sequence``, ``Pair``, and ``Lineage``. These models are heirarchical -- a ``Lineage`` is composed of one 
or more ``Pair`` objects, a ``Pair`` is composed of one or more ``Sequence`` objects -- and contain methods 
appropriate for each level of granularity. 

  * ``Sequence``: model for representing a single antibody sequnce (either heavy or light chain).
    Provides a means to store and access abstar annotations. Includes common methods of sequence
    manipulation, including slicing, reverse-complement, and conversion to FASTA format. The ``Sequence``
    object is used extensively throughout the ab[x] toolkit.
  * ``Pair``: model for representing paired (heavy and light) antibody sequences. Comprised of one 
    or more ``Sequence`` objects. 
  * ``Lineage``: model for representing an antibody clonal lineage. Comprised of one or more ``Pair``
    objects. Includes methods for lineage manipulation, including generating dot alignments and UCA calculation.


tools
------

In addition to the core models, abutils provides a number of commonly used functions. 
These functions are widely used throughput the ab[x] toolkit and are suitable for incorporation 
into custom pipelines or for use when performing interactive analyses:

  * **alignment**: local (Smith-Waterman), global (Needleman-Wunsch) and semi-global pairwise sequence alignment, 
    as well as multiple sequence alignment using MAFFT_ or MUSCLE_

  * **clustering**: identity-based sequence clustering with VSEARCH_, CDHIT_, or MMseqs2_

  * **phylogeny**: computing lineage phylogenies with FastTree_ or IgPhyML_, tree drawing with baltic_



plots
------

``abutils`` provides a




.. _abstar: https://github.com/briney/abstar
.. _MAFFT: https://mafft.cbrc.jp/alignment/software/
.. _MUSCLE: https://www.drive5.com/muscle/
.. _VSEARCH: https://github.com/torognes/vsearch
.. _CDHIT: http://weizhongli-lab.org/cd-hit/
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
.. _FastTree: http://www.microbesonline.org/fasttree/
.. _IgPhyML: https://github.com/kbhoehn/IgPhyML
.. _baltic: https://github.com/evogytis/baltic
