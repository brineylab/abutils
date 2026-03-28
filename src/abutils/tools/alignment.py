#!/usr/bin/env python
# filename: alignment.py
#
# Re-export module: preserves the ``abutils.tools.alignment`` public API
# while the implementation lives in ``msa.py`` and ``pairwise.py``.

from .msa import (  # noqa: F401
    MultipleSequenceAlignment,
    dot_alignment,
    famsa,
    mafft,
    make_consensus,
    muscle,
    muscle_v3,
)
from .pairwise import (  # noqa: F401
    CIGAR,
    CIGARElement,
    GlobalAlignment,
    LocalAlignment,
    PairwiseAlignment,
    SemiGlobalAlignment,
    global_alignment,
    local_alignment,
    process_targets,
    semiglobal_alignment,
)

__all__ = [
    "famsa",
    "mafft",
    "muscle",
    "muscle_v3",
    "local_alignment",
    "global_alignment",
    "semiglobal_alignment",
    "dot_alignment",
    "make_consensus",
    "MultipleSequenceAlignment",
    "PairwiseAlignment",
    "LocalAlignment",
    "GlobalAlignment",
    "SemiGlobalAlignment",
    "CIGAR",
]
