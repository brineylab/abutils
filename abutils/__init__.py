import os

from . import bin, cl, io, log, pl, tl
from . import cl as color
from .core import lineage, pair, sequence
from .core.lineage import Lineage
from .core.pair import Pair
from .core.sequence import Sequence
from .version import __version__

BINARY_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "binaries"))
