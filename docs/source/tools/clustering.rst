.. _clustering:

clustering
=====================

The primary clustering function is :ref:`abutils.tl.cluster <cluster-function>`, which can cluster sequences
using CD-HIT_, VSEARCH_, or MMseqs2_. The function returns a :ref:`abutils.tl.Clusters <clusters-class>` 
object, which contains the clustering results as one or more :ref:`abutils.tl.Cluster <cluster-class>` objects.

examples
-----------


.. code-block:: python

    from abutils.tools import cluster

    # cluster sequences using CD-HIT
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9)

    # cluster sequences using VSEARCH
    clusters = cluster('sequences.fasta', method='vsearch', identity=0.9)

    # cluster sequences using MMseqs2
    clusters = cluster('sequences.fasta', method='mmseqs2', identity=0.9)

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit', identity=0.9, output='clusters.fasta')

    # cluster sequences using CD-HIT, and write the results to a file
    clusters = cluster('sequences.fasta', method='cdhit',


.. _cluster-function:
.. autofunction:: abutils.tools.cluster.cluster

.. _clusters-class:
.. autoclass:: abutils.tools.cluster.Clusters

.. _cluster-class:
.. autoclass:: abutils.tools.cluster.Cluster


.. _CD-HIT: http://weizhongli-lab.org/cd-hit/
.. _VSEARCH: https://github.com/torognes/vsearch
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
