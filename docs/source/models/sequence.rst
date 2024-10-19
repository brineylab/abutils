.. _sequence-class:

sequence
===================

``Sequence`` objects are the fundamental building block for much of the ``abutils`` package. 
Virtually all functions and methods that operate on a one or more sequences will accept ``Sequence`` objects as input.
``Sequence`` can be created from a variety of inputs, including strings, lists, dictionaries, and BioPython_ ``SeqRecord`` 
objects. Below are some brief examples of how to create and use ``Sequence`` objects.

|  

examples
--------

|  


instantiation
................






api
--------


.. autoclass:: abutils.core.sequence.Sequence
    :members:

| 

.. autofunction:: abutils.core.sequence.reverse_complement

| 

.. autofunction:: abutils.core.sequence.translate


.. _biopython: https://biopython.org/wiki/SeqRecord