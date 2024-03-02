.. _clustering:

sequence clustering
===================

The primary clustering function is :ref:`abutils.tl.cluster <cluster-function>`, which can cluster sequences
using CD-HIT_, VSEARCH_, or MMseqs2_. The function returns a :ref:`abutils.tl.Clusters <clusters-class>` 
object, which contains the clustering results as one or more :ref:`abutils.tl.Cluster <cluster-class>` objects.

.. _cluster-function:
.. autofunction:: abutils.tools.cluster.cluster

.. _clusters-class:
.. autoclass:: abutils.tools.cluster.Clusters

.. _cluster-class:
.. autoclass:: abutils.tools.cluster.Cluster


.. _CD-HIT: http://weizhongli-lab.org/cd-hit/
.. _VSEARCH: https://github.com/torognes/vsearch
.. _MMseqs2: https://github.com/soedinglab/MMseqs2
