import os
from importlib.metadata import version

from . import bin, cl, io, log, pl, tl
from . import cl as color
from .core import lineage, pair, sequence
from .core.lineage import Lineage
from .core.pair import Pair
from .core.sequence import Sequence

__version__ = version("abutils")

BINARY_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "binaries"))
