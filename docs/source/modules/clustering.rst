.. _clustering:

sequence clustering
===================

The primary clustering function is :func:`abutils.tl.cluster`, which can cluster sequences
using CD-HIT, VSEARCH, or MMseqs2_. The function returns a :class:`abutils.tl.Clusters` 
object, which contains the clustering results.


.. autofunction:: abutils.tools.cluster.cluster


.. autoclass:: abutils.tools.cluster.Clusters
.. autoclass:: abutils.tools.cluster.Cluster