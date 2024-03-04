.. _pairwise-alignment:

pairwise alignment
===================

``abutils`` can perform several types of pairwise alignments, including global, local, and semi-global alignments.
All pairwise alignment functions return a subclass of the :class:`abutils.tl.PairwiseAlignment` object. The return 
classes for each pairwise alignment type are identical, with the exception of the ``alignment_function`` and ``alignment_type``
properties, so the user can easily switch between alignment methods with minimal changes to their code.


examples
-------------

**Local sequence alignment**
The ``query`` and ``target`` sequences for all ``abutils`` pairwise alignment functions can be an :class:`abutils.Sequence`` 
object, or anything accepted by :class:`abutils.Sequence`. The following example uses sequence strings, 
performs local alignment, and prints the alignment to the console.

.. note::

    If sequences are provided as strings, a random sequence ID will be assigned to each sequence.


.. code-block:: python

    import abutils
    
    # input sequences, as a FASTA-formatted string
    fasta_string = '>seq1\nATGCATGCATGC\n>seq2\nATGCATGCATGC'
    
    # calculate and print the alignment
    aln = abutils.tl.local_alignment(fasta_string)
    print(aln)

| 

**Global sequence alignment with custom scoring parameters**
All ``abutils`` alignment functions accept custom scoring parameters. These parameters are:

  - ``match``: the score for a match (default is ``3``)
  - ``mismatch``: the score for a mismatch (default is ``-2``)
  - ``gap_open``: the penalty for opening a gap (default is ``5``)
  - ``gap_extend``: the penalty for extending a gap (default is ``2``)

To perform a simple identity alignment, set ``match`` to ``1`` and ``mismatch`` to ``0``, like so:

.. code-block:: python

    import abutils
    
    # input sequences, as a FASTA-formatted string
    fasta_string = '>seq1\nATGCATGCATGC\n>seq2\nATGCATGCATGC'
    
    # calculate and print the alignment
    aln = abutils.tl.global_alignment(fasta_string, match=1, mismatch=0)

| 

**Semi-global alignment using a pre-defined similarity matrix**
All ``abutils`` alignment functions can accept any ``parasail`` similarity matrix. The following example uses the ``blosum62`` matrix:

.. code-block:: python

    import abutils
     
    aln = abutils.tl.semi_global_alignment("path/to/sequences.fasta", matrix='blosum62')




| 

api
-----------

.. autofunction:: abutils.tools.alignment.global_alignment
.. autofunction:: abutils.tools.alignment.local_alignment
.. autofunction:: abutils.tools.alignment.dot_alignment

| 

.. autoclass:: abutils.tools.alignment.PairwiseAlignment

| 

.. autoclass:: abutils.tools.alignment.LocalAlignment
.. autoclass:: abutils.tools.alignment.GlobalAlignment
.. autoclass:: abutils.tools.alignment.SemiGlobalAlignment
