![](https://img.shields.io/pypi/v/abutils.svg?colorB=blue)
[![Build Status](https://travis-ci.com/briney/abutils.svg?branch=master)](https://app.travis-ci.com/github/briney/abutils)
[![Documentation Status](https://readthedocs.org/projects/abutils/badge/?version=latest)](https://abutils.readthedocs.io/en/latest/?badge=latest)
![](https://img.shields.io/pypi/pyversions/abutils.svg)
![](https://img.shields.io/badge/license-MIT-blue.svg)

# abutils

Models and general-purpose utilities for working with antibody repertoire data.
abutils is a core component of the ab\[x\] toolkit for antibody sequence analysis.
  
  - Source code: [github.com/briney/abutils](https://github.com/briney/abutils)  
  - Documentation: [abutils.readthedocs.org](http://abutils.readthedocs.org)  
  - Download: [pypi.python.org/pypi/abutils](https://pypi.python.org/pypi/abutils)  
  - Docker: [hub.docker.com/r/briney/abstar/](https://hub.docker.com/r/briney/abstar/)  
  
### install  
`pip install abutils`  


### api  
The intended use of abutils is through the public API, enabling incorporation of abutils' methods and utilities into integrated analysis pipelines, other standalone software tools, or for interative analysis of antibody repertoires. See the abutils [documentation](http://abutils.readthedocs.org) for more detail about the API.  


### testing  
To run the test suite, clone or download the repository and run `pytest ./` from the top-level directory. The same tests are run after every commit using TravisCI.  
  

### requirements  
Python 3.6+  
biopython  
celery  
ete3  
matplotlib  
numpy  
nwalign3  
pandas  
paramiko  
pymongo  
pytest  
scikit bio  
seaborn   
  
All of the above dependencies can be installed with pip, and will be installed automatically when installing abstar with pip.  
If you're new to Python, a great way to get started is to install the [Anaconda Python distribution](https://www.continuum.io/downloads), which includes pip as well as a ton of useful scientific Python packages.
  
abutils has a few additional non-python dependencies that are not required for installation
but are necessary for specific functions:

* ``abutils.alignment.mafft`` requires [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* ``abutils.mongodb.mongoimport`` requires [MongoDB](https://www.mongodb.org/)
* ``abutils.phylogeny.fasttree`` requires [FastTree](http://www.microbesonline.org/fasttree/)
* ``abutils.phylogeny.igphyml`` requires [IgPhyML](https://github.com/kbhoehn/IgPhyML)
* ``abutils.phylogeny.lsd`` [requires LSD](https://github.com/tothuhien/lsd-0.3beta)
* ``abutils.s3`` requires [s3cmd](https://s3tools.org/s3cmd)
