install
=======

The easiest way to install ``abutils`` locally (on macOS or Linux) is to use ``pip``::

    $ pip install abutils

If you don't have ``pip``, the Anaconda_ Python distribution contains ``pip`` along 
with a ton of useful scientific Python packages and is a great way to get 
started with Python.

``abutils`` does not run natively on Windows, but Windows users can use Docker_ 
(the brineylab datascience_ Docker container contains the entire ab[x] toolkit,
which includes ``abutils``)::

    $ docker pull brineylab/datascience
    $ docker run -it brineylab/datascience

Stable_ and development_ versions of ``abutils`` can also be downloaded from Github. 
You can manually install the latest version of ``abutils`` with::

    $ git clone https://github.com/briney/abutils
    $ cd abutils/
    $ python setup.py install


requirements
------------

* Python 3.7+
* baltic_
* biopython_
* celery_
* ete3_
* fastcluster_
* matplotlib_
* mnemonic_
* numpy_
* nwalign3_
* pandas_
* paramiko_
* parasail_
* pymongo_
* pytest_
* python-circos_
* python-Levenshtein_
* pyyaml_
* sample-sheet_
* scikit-learn_
* scipy_
* seaborn_


additional dependencies
-----------------------

Whenever possible, ``abutils`` bundles required third-party binaries, but there are a few 
additional non-python dependencies that must be separately installed. These tools are 
not needed for installation, but are necessary for specific functions:

* ``abutils.tl.igphyml`` requires IgPhyML_
* ``abutils.tl.lsd`` requires LSD_
* ``abutils.s3`` requires s3cmd_

If using Docker, all of the the optional dependencies are included.


.. _Docker: https://www.docker.com/
.. _datascience: https://hub.docker.com/repository/docker/brineylab/datascience/general
.. _Anaconda: https://www.continuum.io/downloads
.. _stable: https://github.com/briney/abstar/releases
.. _development: https://github.com/briney/abstar
.. _abutils: https://github.com/briney/abutils
.. _biopython: http://biopython.org/
.. _celery: http://www.celeryproject.org/
.. _scikit bio: http://scikit-bio.org/
.. _pymongo: https://api.mongodb.org/python/current/
.. _MongoDB: https://www.mongodb.org/
.. _pytest: https://docs.pytest.org/en/latest/
.. _ete3: http://etetoolkit.org/
.. _matplotlib: https://matplotlib.org/
.. _numpy: http://www.numpy.org/
.. _nwalign3: https://github.com/briney/nwalign3
.. _pandas: https://pandas.pydata.org/
.. _paramiko: http://www.paramiko.org/
.. _seaborn: https://seaborn.pydata.org/
.. _MAFFT: https://mafft.cbrc.jp/alignment/software/
.. _s3cmd: https://s3tools.org/s3cmd
.. _FastTree: http://www.microbesonline.org/fasttree/
.. _IgPhyML: https://github.com/kbhoehn/IgPhyML
.. _LSD: https://github.com/tothuhien/lsd-0.3beta
.. _parasail: https://github.com/jeffdaily/parasail-python
.. _baltic: https://github.com/evogytis/baltic
.. _fastcluster: https://github.com/dmuellner/fastcluster
.. _mnemonic: https://github.com/trezor/python-mnemonic
.. _python-circos: https://github.com/ponnhide/pyCircos
.. _python-Levenshtein: https://github.com/ztane/python-Levenshtein
.. _pyyaml: https://pyyaml.org/
.. _sample-sheet: https://github.com/clintval/sample-sheet
.. _scikit-learn: https://scikit-learn.org/stable/
.. _scipy: https://www.scipy.org/
