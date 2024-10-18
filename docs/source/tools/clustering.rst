.. _clustering:

clustering
=====================

The primary clustering function is ``abutils.tl.cluster``, which can cluster sequences
using CD-HIT_, VSEARCH_, or MMseqs2_. The function returns a :class:`abutils.tl.Clusters` 
object, which contains the clustering results as one or more :class:`abutils.tl.Cluster` objects.

``abutils.tl.cluster`` can accept a variety of inputs, including:

- a path to a FASTA file
- a FASTA-formatted string
- a list of ``abutils.Sequence`` objects
- a list of anything accepted by :class:`abutils.Sequence`

The ``threshold`` argument is the sequence identity threshold for clustering, and should be between 0.0 and 1.0.

The ``algo`` argument selects the clustering algorithm to use. It can be ``'cdhit'``, ``'vsearch'``, or ``'mmseqs2'``. 
If ``algo`` is not provided, ``abutils.tl.cluster`` will use CD-HIT for inputs with fewer than 1000 sequences or MMseqs2 for inputs with 1000 or more sequences.
Binaries for each of the available clustering algorithms are packaged with ``abutils``, however, if you would like to use an alternate binary (for example, to use 
a different version of one of the built-in binaries), you can specify the path to the desiredbinary using the following optional arguments:

- ``cdhit_bin``
- ``vsearch_bin``
- ``mmseqs_bin``

|  

examples
-----------

|  

**clustering with CD-HIT at 90% identity, using a FASTA file as input**  
  
The ``sequences`` argument can be a path to a FASTA file, a FASTA-formatted string,
a list of ``Sequence`` objects, or a list of anything accepted by 
:class:`abutils.Sequence`. The ``algo`` argument can be ``'cdhit'``, ``'vsearch'``, or ``'mmseqs2'``. 
If ``algo`` is not provided, ``cluster()`` will use CD-HIT for inputs with fewer than 1000 
sequences or MMseqs2 for inputs with 1000 or more sequences. The ``threshold`` argument 
is the sequence identity threshold for clustering.
  
.. code-block:: python

    import abutils

    clusters = abutils.tl.cluster(
        sequences='path/to/sequences.fasta', 
        algo='cdhit', 
        threshold=0.9
    )

| 

**get the largest cluster from VSEARCH clustering of a list of Sequence objects**
  
Iterating over a :class:`abutils.tl.Clusters` object iterates over all of the clusters it contains, 
which are themselves :class:`abutils.tl.Cluster` objects. By default, they are sorted by size in 
descending order, so the first cluster is the largest. 
  
.. code-block:: python

    import abutils

    sequences = abutils.io.read_fasta('path/to/sequences.fasta')
    clusters = abutils.tl.cluster(
        sequences=sequences, 
        algo='vsearch', 
        threshold=0.9
    )
    largest_cluster = clusters[0]
  
| 

**calculate the consensus sequence of the largest cluster, using MMSeqs2**
  
:class:`abutils.tl.Cluster` objects have a ``consensus`` property that returns an :class:`abutils.Sequence` object
representing the consensus sequence of the cluster. The ``consensus`` property is lazy, 
meaning it is not calculated until it is accessed, and once calculated, it is cached. 
Under the hood, the consensus sequence is calculated using the ``make_consensus()`` 
method, which can also be called directly to provides more control over the consensus 
sequence generation process. Calling ``make_consensus()`` automatically saves the consensus 
sequence to the ``consensus`` property (and overwrites any cached consensus sequence).

.. code-block:: python

    import abutils

    clusters = abutils.tl.cluster(
        sequences='path/to/sequences.fasta', 
        algo='mmseqs2', 
        threshold=0.9
    )
    largest_cluster = clusters[0]
    consensus = largest_cluster.consensus

    # to overwrite the cached consensus sequence
    largest_cluster.make_consensus()

| 

api
-----------

.. _cluster-function:
  
.. autofunction:: abutils.tl.cluster

| 

.. _clusters-class:
  
.. autoclass:: abutils.tl.Clusters
    :members:
| 

.. _cluster-class:
  
.. autoclass:: abutils.tl.Cluster
    :members:

.. _CD-HIT: http://weizhongli-lab.org/cd-hit/
.. _VSEARCH: https://github.com/torognes/vsearch
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
