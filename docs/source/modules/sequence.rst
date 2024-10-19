.. _sequence-class:

sequence
===================

``Sequence`` objects are the fundamental building block for much of the ``abutils`` package. 
Virtually all functions and methods that operate on a one or more sequences will accept ``Sequence`` objects as input.
``Sequence`` can be created from a variety of inputs, including strings, lists, dictionaries, and BioPython_ ``SeqRecord`` 
objects. Below are some brief examples of how to create and use ``Sequence`` objects.

|  

instantiation
----------------

``abutils`` has a number of convenience functions for batch creation of ``Sequence`` objects 
common file formats, including FASTA, FASTQ, AIRR-C_, and Parquet. Details and examples of these 
functions can be found in the :ref:`io` section.

Individual ``Sequence`` objects can be created from a string:

.. code-block:: python

    import abutils

    # create a sequence from a string
    sequence = abutils.Sequence("ATCG")


 .. note::

    If provided a string, the sequence ID will be randomly generated if not specified.
    To specify the sequence ID, you can pass a ``Sequence`` object to the ``id`` argument:

    .. code-block:: python

        # create a sequence from a string
        sequence = abutils.Sequence("ATCG", id="my_sequence")

|  

``Sequence`` objects can also be created from a list, of the form ``[id, sequence]``:

.. code-block:: python

    # create a sequence from a list
    sequence = abutils.Sequence(["my_sequence", "ATCG"])

|  

``Sequence`` objects can also be created from a dictionary, which provides a means for 
including additoinal annotations beyond just the sequence and ID:

.. code-block:: python

    # create a sequence from a dictionary
    sequence = abutils.Sequence({"sequence_id": "my_sequence", "sequence": "ATCG", "productive": True})

    # all annotations can be accessed using dictionary-style indexing
    sequence["productive"]

.. note::

    Dictionary keys are typically expected to follow the naming conventions of the 
    AIRR-C_ rearrangement schema. The ``Sequence`` object will automatically populate the special 
    properties ``id`` and ``sequence`` from the provided dictionary if the correct key names (``"sequence_id"`` 
    and ``"sequence"``, respectively) are used.

|  

usage
----------------

``Sequence`` objects have several convenient properties for common sequence manipulations:

.. code-block:: python

    # reverse complement
    rc = sequence.reverse_complement()

    # translate
    aa = sequence.translate()

|  





api
--------


.. autoclass:: abutils.core.sequence.Sequence
    :members:

| 

.. autofunction:: abutils.core.sequence.reverse_complement

| 

.. autofunction:: abutils.core.sequence.translate


.. _biopython: https://biopython.org/wiki/SeqRecord
.. _AIRR-C: https://docs.airr-community.org/en/latest/datarep/rearrangements.html