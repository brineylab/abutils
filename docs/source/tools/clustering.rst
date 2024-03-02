.. _clustering:

clustering
=====================

The primary clustering function is :ref:`abutils.tl.cluster <cluster-function>`, which can cluster sequences
using CD-HIT_, VSEARCH_, or MMseqs2_. The function returns a :ref:`abutils.tl.Clusters <clusters-class>` 
object, which contains the clustering results as one or more :ref:`abutils.tl.Cluster <cluster-class>` objects.

examples
-----------

**clustering with CD-HIT at 90% identity, using a FASTA file as input**  
  
The ``sequences`` argument can be a path to a FASTA file, a FASTA-formatted string,
a list of ``Sequence`` objects, or a list of :ref:`anything accepted by <sequence-class>` 
``Sequence``. The ``algo`` argument can be ``'cdhit'``, ``'vsearch'``, or ``'mmseqs2'``. 
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


**get the largest cluster from VSEARCH clustering of a list of ``Sequence`` objects**
  
Iterating over a ``Clusters`` object iterates over all of the clusters it contains, 
which are themselves ``Cluster`` objects. By default, they are sorted by size in 
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





.. _cluster-function:
.. autofunction:: abutils.tools.cluster.cluster

.. _clusters-class:
.. autoclass:: abutils.tools.cluster.Clusters

.. _cluster-class:
.. autoclass:: abutils.tools.cluster.Cluster


.. _CD-HIT: http://weizhongli-lab.org/cd-hit/
.. _VSEARCH: https://github.com/torognes/vsearch
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
