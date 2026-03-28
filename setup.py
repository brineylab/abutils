import os

from setuptools import find_packages, setup

# read version
version_file = os.path.join(os.path.dirname(__file__), "abutils", "version.py")
with open(version_file) as f:
    exec(f.read())

# read requirements
with open("requirements.txt") as f:
    requirements = f.read().splitlines()

# read long description
with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="abutils",
    version=__version__,
    author="Bryan Briney",
    author_email="briney@scripps.edu",
    description="Utilities for analysis of adaptive immune receptor repertoire (AIRR) data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/briney/abutils",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.10",
)
