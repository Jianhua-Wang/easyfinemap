"""Perform Fine-mapping.

Perform fine-mapping for a locus using the following methods:
1. LD-free
2. LD-based
    2.1. without annotation
        2.1.1. FINEMAP
        2.1.2. PAINTOR
        2.1.3. CAVIARBF
        2.1.4. SuSiE
    2.2. with annotation
        2.2.1. PolyFun+SuSiE
"""

import logging
from pathlib import Path

from easyfinemap.ldref import LDRef
from easyfinemap.loci import Loci
from easyfinemap.sumstat import SumStat
from easyfinemap.tools import Tools
from easyfinemap.utils import get_significant_snps, io_in_tempdir, make_SNPID_unique


class EasyFinemap(object):
    """Main class."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger('EasyFinemap')
        tool = Tools()
        self.finemap = tool.finemap
        self.paintor = tool.paintor
        self.gcta = tool.gcta
        self.plink = tool.plink
        self.bcftools = tool.bcftools
        self.caviarbf = tool.caviarbf
        self.tmp_root = Path.cwd() / "tmp" / "finemapping"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)

    def run(self):
        """Run the program."""
        self.logger.info("Running the program")

    @io_in_tempdir('./tmp/finemapping')
    def make_ld(self, sumstat: SumStat, ldref: LDRef, out: str):
        """Make the LD matrix."""
        self.logger.info("Making the LD matrix")
