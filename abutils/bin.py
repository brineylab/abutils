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

__all__ = [
    "cdhit",
    "fastp",
    "fasttree",
    "mmseqs",
    "mafft",
    "muscle",
    "muscle_v3",
    "vsearch",
]


BIN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bin")
SYSTEM = platform.system().lower()
MACHINE = platform.machine().lower().replace("x86_64", "amd64")

cdhit = os.path.join(BIN_DIR, f"cd-hit_{SYSTEM}_{MACHINE}")
fastp = os.path.join(BIN_DIR, f"fastp_{SYSTEM}")
fasttree = os.path.join(BIN_DIR, f"fasttree_{SYSTEM}_{MACHINE}")
mafft = os.path.join(BIN_DIR, f"mafft_{SYSTEM}_{MACHINE}/mafft.bat")
mmseqs = os.path.join(BIN_DIR, f"mmseqs_{SYSTEM}_{MACHINE}")
muscle = os.path.join(BIN_DIR, f"muscle_{SYSTEM}_{MACHINE}")
muscle_v3 = os.path.join(BIN_DIR, f"muscle3_{SYSTEM}_{MACHINE}")
vsearch = os.path.join(BIN_DIR, f"vsearch_{SYSTEM}_{MACHINE}")
