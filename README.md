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
  - Docker: [hub.docker.com/r/brineylab/datascience/](https://hub.docker.com/r/brineylab/datascience/)  
  
## install  
``` bash
pip install abutils
```


## api  
We've tried to design the  `abutils` API to be intuitive yet powerful, with the goal of enabling both interactive analyses (via environments like Jupyter notebooks) as well as integration of `abutils` tools into more complex analysis pipelines and/or standalone software tools. See the [documentation](http://abutils.readthedocs.org) for more detail about the API. As always, any feedback is greatly appreciated!!  


### testing  
You can run the complete `abutils` test suite by first installing `pytest`:
``` bash
pip install pytest
```

and then running:

``` bash
git clone https://github.com/brineylab/abutils
cd abutils
pytest
```

This test suite is automatically run against all supported versions of Python following every commit.
  

### requirements  
**python 3.10+**  
  
abstar  
baltic  
biopython  
dnachisel  
fastcluster  
matplotlib  
mnemonic  
natsort  
numpy  
pandas  
parasail  
polars  
prettytable  
pyarrow  
pyfamsa  
pyfastx  
pytest  
python-circos  
pyyaml  
rapidfuzz  
sample-sheet  
scikit-learn  
scipy  
seaborn  
smart_open  
tqdm  
  
`abutils` includes several additional binaries that are required for certain functionality:

* ``abutils.tl.mafft`` uses [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* ``abutils.tl.muscle`` uses [MUSCLE](https://www.drive5.com/muscle/)
* ``abutils.tl.cluster`` uses [CD-HIT](https://cd-hit.org), [MMseqs2](https://github.com/soedinglab/MMseqs2), and [VSEARCH](https://github.com/torognes/vsearch)
* ``abutils.tl.fasttree`` uses [FastTree](http://www.microbesonline.org/fasttree/)

Although these binaries are all packaged into `abutils`, each respective `abutils.tl` function provides the option to supply a alternate binary path in case you'd prefer to use a different version.  


