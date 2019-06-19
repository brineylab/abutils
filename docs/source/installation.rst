install
=======

The easiest way to install abutils locally (on macOS or Linux) is to use pip::

    $ pip install abutils

If you don't have pip, the Anaconda_ Python distribution contains pip along 
with a ton of useful scientific Python packages and is a great way to get 
started with Python.

abutils does not run natively on Windows, but Windows users can use Docker_ (the abstar Docker container also includes abutils)::

    $ docker pull briney/abstar
    $ docker run -it briney/abstar

Stable_ and development_ versions of abstar can also be downloaded from Github. 
You can manually install the latest version of abstar with::

    $ git clone https://github.com/briney/abutils
    $ cd abutils/
    $ python setup.py install

.. note::

    If installing manually via setup.py and you don't already have scikit-bio installed, 
    you may get an error when setuptools attempts to install scikit-bio. This can be fixed 
    by first installing scikit-bio with pip::

        $ pip install scikit-bio

    and then retrying the manual install of abutils. Starting with version 0.5, scikit-bio 
    dropped support for Python 2.7, so install scikit-bio on Python 2.7 with::

        $ pip install scikit-bio<=0.4.2


requirements
------------

* Python 3.5+
* biopython_
* celery_
* ete3_
* matplotlib_
* numpy_
* nwalign3_
* pandas_
* paramiko_
* pymongo_
* pytest_
* `scikit bio`_
* seaborn_


additional dependencies
---------------------

abutils has a few additional non-python dependencies that are not required for installation
but are necessary for specific functions:

* ``abutils.alignment.mafft`` requires MAFFT_
* ``abutils.mongodb.mongoimport`` requires MongoDB_
* ``abutils.phylogeny.fasttree`` requires FastTree_
* ``abutils.phylogeny.igphyml`` requires IgPhyML_
* ``abutils.phylogeny.lsd`` requires LSD_
* ``abutils.s3`` requires s3cmd_

If using Docker, all of the the optional dependencies are included.


.. _Docker: https://www.docker.com/
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




