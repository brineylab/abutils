#!/usr/bin/env python
# filename: bin.py
#
#
# Copyright (c) 2023 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


import os
import platform
from typing import Optional
from urllib.request import urlretrieve

__all__ = ["get_path", "get_binary_directory"]

# __all__ = [
#     "cdhit",
#     "fastp",
#     "fasttree",
#     "mmseqs",
#     "mafft",
#     "muscle",
#     "muscle_v3",
#     "vsearch",
# ]


BIN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binaries")
SYSTEM = platform.system().lower()
MACHINE = platform.machine().lower().replace("x86_64", "amd64")

# cdhit = os.path.join(BIN_DIR, f"cd-hit_{SYSTEM}_{MACHINE}")
# fastp = os.path.join(BIN_DIR, f"fastp_{SYSTEM}")
# fasttree = os.path.join(BIN_DIR, f"fasttree_{SYSTEM}_{MACHINE}")
# mafft = os.path.join(BIN_DIR, f"mafft_{SYSTEM}_{MACHINE}/mafft.bat")
# mmseqs = os.path.join(BIN_DIR, f"mmseqs_{SYSTEM}_{MACHINE}")
# muscle = os.path.join(BIN_DIR, f"muscle_{SYSTEM}_{MACHINE}")
# muscle_v3 = os.path.join(BIN_DIR, f"muscle3_{SYSTEM}_{MACHINE}")
# vsearch = os.path.join(BIN_DIR, f"vsearch_{SYSTEM}_{MACHINE}")


def get_binary_directory() -> str:
    """
    Gets the abutils binary directory.

    Returns
    -------
    str
        Path to the abutils binary directory.

    """
    return BIN_DIR


def copy_to_binary_directory(source: str, name: Optional[str] = None) -> None:
    """
    Copies a binary to the abutils binary directory.

    Parameters
    ----------
    source : str
        Path to the binary to be copied.

    name : str, optional
        Name of the binary in the binary directory. If not provided, the basename of the
        source binary will be used.

    """
    binary_dir = get_binary_directory()
    name = name if name is not None else os.path.basename(source)
    target = os.path.join(binary_dir, source)
    if not os.path.exists(target):
        os.system(f'cp "{source}" "{target}"')


def get_path(binary: str) -> str:
    """
    Gets the path to a pre-packaged binary.

    .. note::
        If the binary is not found in the binary directory, it will be downloaded from
        the brineylab S3 bucket. This is because some binaries are too large
        to be included in the ``abutils`` PyPI package.

    Parameters
    ----------
    binary : str
        Name of the binary. Available binaries are:
            - blastn
            - cd-hit
            - fastp
            - fasttree
            - makeblastdb
            - mafft
            - minimap2
            - mmseqs
            - muscle
            - muscle_v3
            - sickle
            - vsearch

    Returns
    -------
    str
        Path to the binary.

    """
    #  format (and maybe fix) binary name
    binary = binary.lower()
    if binary == "mmseqs2":
        binary = "mmseqs"
    if binary == "muscle_v3":
        binary = "muscle3"
    if binary == "cd-hit":
        binary = "cdhit"
    if binary == "minimap":
        binary = "minimap2"
    available = [
        # blast
        "blastn",
        "makeblastdb",
        # clustering/alignment
        "cdhit",
        "mafft",
        "minimap2",
        "mmseqs",
        "muscle",
        "muscle3",
        "vsearch",
        # qc/preprocessing
        "fastp",
        "sickle",
        # phylogeny
        "fasttree",
    ]
    if binary not in available:
        raise ValueError(
            f"Binary '{binary}' not available. Available binaries: \n"
            + "\n".join(available)
        )
    # get binary path
    if binary in ["fastp", "muscle3", "blastn", "makeblastdb"]:
        bin_name = f"{binary}_{SYSTEM}"
        bin_path = os.path.join(BIN_DIR, bin_name)
    elif binary == "mafft":
        bin_name = f"{binary}_{SYSTEM}/mafft.bat"
        bin_path = os.path.join(BIN_DIR, bin_name)
    else:
        bin_name = f"{binary}_{SYSTEM}_{MACHINE}"
        bin_path = os.path.join(BIN_DIR, bin_name)
    # download binary if missing
    if not os.path.exists(bin_path):
        print(
            f"Downloading missing binary: {binary}. This will only occur once on intial use of the binary."
        )
        download_missing_binary(bin_name, bin_path)
    return bin_path


# download missing binary (some are too large to included in the PyPI package)
def download_missing_binary(bin_name: str, bin_path: str) -> None:
    url = f"https://brineylab.s3.amazonaws.com/tools/abutils_binaries/{bin_name}"
    urlretrieve(url, bin_path)


# # download any missing binaries (some are too large to included in the PyPI package)
# def download_missing_binaries():
#     to_download = []
#     for b in [
#         "cd-hit",
#         "fastp",
#         "fasttree",
#         "mafft",
#         "mmseqs",
#         "muscle",
#         "muscle3",
#         "vsearch",
#     ]:
#         if b in ["fastp", "muscle3"]:
#             bin_name = f"{b}_{SYSTEM}"
#             bin_path = os.path.join(BIN_DIR, bin_name)
#         else:
#             bin_name = f"{b}_{SYSTEM}_{MACHINE}"
#             bin_path = os.path.join(BIN_DIR, bin_name)
#         if not os.path.exists(bin_path):
#             to_download.append((bin_name, bin_path))
#     if to_download:
#         print(
#             "Downloading missing binaries. This will only occur once on intial use of abutils."
#         )
#         url = f"https://brineylab.s3.amazonaws.com/tools/abutils_binaries/{bin_name}"
#         urlretrieve(url, bin_path)


# download_missing_binaries()
