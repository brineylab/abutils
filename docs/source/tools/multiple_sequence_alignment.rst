.. _multiple-sequence-alignment:

multiple sequence alignment
===========================

``abutils`` can perform multiple sequence alignments using either MAFFT_, MUSCLE_, or FAMSA_ (using the pyfamsa_ package). 
All multiple sequence alignment functions return a :class:`abutils.tl.MultipleSequenceAlignment` object, which builds on 
the :class:`Bio.Align.MultipleSeqAlignment` class from Biopython. This object provides a number of
convenient methods for working with the alignment, including writing to file, trimming, and calculating
consensus sequences. Additionally, because the same object is returned regardless of the alignment method used,
the user can easily switch between alignment methods with minimal changes to their code.

|  

.. csv-table:: 
   :header: "tool", "function"
   :widths: 10, 16

   "MAFFT", :ref:`abutils.tl.mafft() <mafft-aln>`
   "MUSCLE", :ref:`abutils.tl.muscle() <muscle-aln>`
   "FAMSA", :ref:`abutils.tl.famsa() <famsa-aln>`


examples
---------

**multiple sequence alignment with MAFFT**

Each of the multiple seqeunce alignment functions can accept a path to a FASTA file, a FASTA-formatted string,
a list of :class:`abutils.Sequence` objects, or a list of anything accepted by :class:`abutils.Sequence`. By default,
calling the alignment function with a list of sequences will return a :class:`abutils.tl.MultipleSequenceAlignment`.


.. code-block:: python

    import abutils

    msa = abutils.tl.mafft('path/to/sequences.fasta')

| 

**multiple sequence alignment with MUSCLE, with the results written to file**

Rather than returning a :class:`abutils.tl.MultipleSequenceAlignment` object, the user can specify a path to which
the alignment file should be written. This is done by passing a file path to the ``alignment_file`` argument.

.. code-block:: python

    import abutils

    # read in a fasta file
    seqs = abutils.io.read_fasta('path/to/sequences.fasta')

    # align the sequences using MUSCLE and write the results to a file
    abutils.tl.muscle(
        seqs, 
        alignment_file='path/to/alignment.fasta', 
        as_file=True
    )

| 

**using a custom binary for multiple sequence alignment**

``abutils`` packages binaries for both MAFFT and MUSCLE, meaning these packages don't need to be separately installed
by the user. However, both :ref:`abutils.tl.mafft() <mafft-aln>` and :ref:`abutils.tl.muscle() <muscle-aln>` allow the 
user to specify the path to a custom binary if desired. For MAFFT, this is done using the ``mafft_bin`` argument, and 
for MUSCLE, the ``muscle_bin`` argument. Both accept a path to the binary as a string.

.. code-block:: python

    import abutils

    # align the sequences using a custom MAFFT binary
    mafft_msa = abutils.tl.mafft(
        'path/to/sequences.fasta', 
        mafft_bin='path/to/mafft'
    )

    # align the sequences using a custom MUSCLE binary
    muscle_msa = abutils.tl.muscle(
        'path/to/sequences.fasta', 
        muscle_bin='path/to/muscle'
    )


| 

api
-------

.. _mafft-aln:

.. autofunction:: abutils.tl.mafft

.. _muscle-aln:

.. autofunction:: abutils.tl.muscle

.. _famsa-aln:

.. autofunction:: abutils.tl.famsa

| 

.. autoclass:: abutils.tl.MultipleSequenceAlignment
    :members:



.. _MAFFT: https://mafft.cbrc.jp/alignment/software/

.. _MUSCLE: https://www.drive5.com/muscle/

.. _FAMSA: https://github.com/refresh-bio/FAMSA

.. _pyfamsa: https://github.com/althonos/pyfamsa
