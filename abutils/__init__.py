# from .core import *
# from .utils import *

import os

from . import cl, io, pl, tl
from . import cl as color
from .core import lineage, pair, sequence
from .core.lineage import Lineage
from .core.pair import Pair
from .core.sequence import Sequence
from .tools.phylo import Phylogeny
from .utils import (
    alignment,
    cluster,
    # color,
    decorators,
    jobs,
    log,
    phylogeny,
    pipeline,
    progbar,
    s3,
    seqio,
    utilities,
)
from .utils import alignment as aln
from .utils import pipeline as path

# from .utils.alignment import SSWAlignment, NWAlignment
from .utils.alignment import GlobalAlignment, LocalAlignment, SemiGlobalAlignment
from .utils.jobs import *
from .utils.progbar import progress_bar
from .version import __version__

BINARY_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "bin"))
