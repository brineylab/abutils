.. _multiple-sequence-alignment:

multiple sequence alignment
===========================

``abutils`` can perform multiple sequence alignments using either MAFFT_ or MUSCLE_. 
All alignment functions return a :class:`abutils.tools.alignment.MultipleSequenceAlignment` object, which builds on 
the :class:`Bio.Align.MultipleSeqAlignment` class from Biopython. This object provides a number of
convenient methods for working with the alignment, including writing to file, trimming, and calculating
consensus sequences. Additionally, because the same object is returned regardless of the alignment method used,
the user can easily switch between alignment methods with minimal changes to their code.


examples
---------




.. autofunction:: abutils.tl.mafft
.. autofunction:: abutils.tools.alignment.muscle

.. autoclass:: abutils.tools.alignment.MultipleSequenceAlignment





.. _MAFFT: https://mafft.cbrc.jp/alignment/software/
.. _MUSCLE: https://www.drive5.com/muscle/
