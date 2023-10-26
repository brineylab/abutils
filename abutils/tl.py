#!/usr/bin/env python
# filename: tl.py


#
# Copyright (c) 2020 Bryan Briney
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


from .core.lineage import group_lineages
from .core.pair import assign_pairs
from .core.sequence import translate

from .tools.alignment import *
from .tools.clonify import clonify
from .tools.cluster import *
from .tools.phylo import *
from .tools.similarity import repertoire_similarity

# from .utils.alignment import (
#     global_alignment,
#     local_alignment,
#     dot_alignment,
#     muscle,
#     muscle_v3,
#     mafft,
# )

# from .utils.cluster import cluster, cluster_vsearch, cluster_mmseqs
from .utils.progbar import progress_bar
from .utils.phylogeny import lsd, igphyml
