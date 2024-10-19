.. _pairwise-alignment:

pairwise alignment
===================

``abutils`` can perform several types of pairwise alignments, including global, local, and semi-global alignments.
All pairwise alignment functions return a subclass of the :class:`abutils.tl.PairwiseAlignment` object. The return 
classes for each pairwise alignment type are identical, with the exception of the ``alignment_function`` and ``alignment_type``
properties, so the user can easily switch between alignment methods with minimal changes to their code.


examples
-------------

**local sequence alignment**  
  
The ``query`` and ``target`` sequences for all ``abutils`` pairwise alignment functions can be an :class:`abutils.Sequence`` 
object, or anything accepted by :class:`abutils.Sequence`. The following example uses sequence strings, 
performs local alignment, and prints the alignment to the console.

.. note::

    If sequences are provided as strings, a random sequence ID will be assigned to each sequence.


.. code-block:: python

    import abutils
    
    # input sequences, as strings
    seq1 = 'ATGCATGCATGC'
    seq2 = 'ATGCATGCATGC'
    
    # calculate and print the alignment
    aln = abutils.tl.local_alignment(seq1, seq2)
    print(aln)

| 

**global sequence alignment with custom scoring parameters**  

All ``abutils`` pairwise alignment functions accept custom scoring parameters. These parameters are:

  - ``match``: the score for a match (default is ``3``)
  - ``mismatch``: the score for a mismatch (default is ``-2``)
  - ``gap_open``: the penalty for opening a gap (default is ``5``)
  - ``gap_extend``: the penalty for extending a gap (default is ``2``)

Alignment functions can also accept any ``parasail`` similarity matrix, passed 
to the ``matrix`` argument. The following example uses the ``blosum62`` matrix:

.. code-block:: python

    import abutils
    
    # input sequences, as strings
    seq1 = 'ATGCATGCATGC'
    seq2 = 'ATGCATGCATGC'
    
    # calculate and print the alignment
    aln = abutils.tl.global_alignment(
        seq1, 
        seq2, 
        matrix='blosum62',
    )

| 

**semi-global alignment against multiple targets and selecting the best match**  

All ``abutils`` pairwise alignment functions can align the same ``query`` sequence against multiple ``target`` sequences
using the ``targets`` argument. This provides a moderate speed increase and avoids the need to loop through the targets. 
The following example aligns a single query sequence against multiple target sequences, sorts the resulting list of 
alignments (which, by default, sorts by alignment score), and selects the top scoring alignment.

.. code-block:: python

    import abutils
    
    # query and target sequences
    query = 'ATGCATGCATGC'
    targets = abutils.io.read_fasta('path/to/targets.fasta')

    alns = abutils.tl.semi_global_alignment(
        query,
        targets=targets
    )

    # get the highest scoring alignment
    best_aln = sorted(alns, reverse=True)[0]


| 

api
-----------

.. autofunction:: abutils.tl.local_alignment
.. autofunction:: abutils.tl.global_alignment
.. autofunction:: abutils.tl.semiglobal_alignment

| 

.. autoclass:: abutils.tl.PairwiseAlignment
    :members:
    :special-members: __init__

| 

.. autoclass:: abutils.tl.LocalAlignment
.. autoclass:: abutils.tl.GlobalAlignment
.. autoclass:: abutils.tl.SemiGlobalAlignment
