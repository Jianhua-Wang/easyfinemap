"""Top-level package for easy_finemap."""

import logging

from rich.logging import RichHandler

from .easyfinemap import EasyFinemap
from .ldref import LDRef
from .loci import Loci
from .plots import locus_plot

# from .sumstat import SumStat

__author__ = """Jianhua Wang"""
__email__ = 'jianhua.mert@gmail.com'
__version__ = '0.4.0'


logging.basicConfig(
    level=logging.WARNING,
    format="%(name)s - %(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)
