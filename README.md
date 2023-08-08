![](https://img.shields.io/pypi/v/abutils.svg?colorB=blue)
[![tests](https://github.com/briney/abutils/actions/workflows/pytest.yml/badge.svg)](https://github.com/briney/abutils/actions/workflows/pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/abutils/badge/?version=latest)](https://abutils.readthedocs.io/en/latest/?badge=latest)
![](https://img.shields.io/pypi/pyversions/abutils.svg)
![](https://img.shields.io/badge/license-MIT-blue.svg)

# abutils

Models, functions and visualization tools for working with adaptive immune receptore repertoire (AIRR) data. The primary purpose of `abutils` is to provide generalizable tools suitable for direct use analyzing bulk AIRR datasets, and is used by [`scab`](https://github.com/briney/scab) for single cell AIRR analysis. `abutils` is a core component of the ab\[x\] toolkit for AIRR data analysis.
  
  - Source code: [github.com/briney/abutils](https://github.com/briney/abutils)  
  - Documentation: [abutils.readthedocs.org](http://abutils.readthedocs.org)  
  - Download: [pypi.python.org/pypi/abutils](https://pypi.python.org/pypi/abutils)  
  <!-- - Docker: [hub.docker.com/r/briney/abstar/](https://hub.docker.com/r/briney/abstar/)   -->  
  
### install  
`pip install abutils`  


### api  
We've tried to design the  `abutils` API to be intuitive yet powerful, with the goal of enabling both interactive analyses (via environments like Jupyter notebooks) as well as integration of `abutils` tools into more complex analysis pipelines and/or standalone software tools. See the [documentation](http://abutils.readthedocs.org) for more detail about the API. As always, any feedback is greatly appreciated!!  


### testing  
You can run the complete `abutils` test suite by first installing `pytest`:
```
pip install pytest
```

followed by:

```
git clone https://github.com/briney/abutils
cd abutils
pytest
```

This test suite is automatically run following every commit, and is tested against all supported versions of Python.
  

### requirements  
**python 3.8+**  
  
abstar  
baltic  
biopython  
celery  
ete3  
fastcluster  
matplotlib  
mnemonic  
natsort  
numpy  
pandas  
paramiko  
parasail  
pytest  
python-circos  
python-Levenshtein  
pyyaml  
sample-sheet  
scikit-learn  
scipy  
seaborn  
smart_open  
  
All of the above dependencies can be installed with `pip`, and will be installed automatically when installing `abutils` with `pip`.  
  
`abutils` packages several additional external binaries that are required for specific functions:

* ``abutils.tl.mafft`` uses [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* ``abutils.tl.muscle`` uses [MUSCLE](https://www.drive5.com/muscle/)
* ``abutils.tl.cluster`` requires:
  * [CD-HIT](https://cd-hit.org)
  * [MMseqs2](https://github.com/soedinglab/MMseqs2)
  * [VSEARCH](https://github.com/torognes/vsearch)
* ``abutils.tl.fasttree`` requires [FastTree](http://www.microbesonline.org/fasttree/)

Althogh these binaries are all packaged into `abutils`, each respective `abutils` function provides the option to supply a different binary path in case you'd prefer to use a different version or an alternate compilation.  


