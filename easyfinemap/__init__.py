"""Top-level package for easy_finemap."""

import logging

from rich.logging import RichHandler

from .easyfinemap import EasyFinemap
from .ldref import LDRef
from .loci import Loci

# from .sumstat import SumStat

__author__ = """Jianhua Wang"""
__email__ = 'jianhua.mert@gmail.com'
__version__ = '0.2.3'


logging.basicConfig(
    level=logging.NOTSET,
    format="%(name)s - %(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)
