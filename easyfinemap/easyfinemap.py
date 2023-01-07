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
from subprocess import PIPE, run
from typing import List, Optional

import pandas as pd

from easyfinemap.constant import ColName
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

    @io_in_tempdir('./tmp/finemapping')
    def cojo_cond(
        self,
        sumstats: pd.DataFrame,
        ldref: str,
        cond_snps: List[str],
        out: str,
        sample_size: int,
        use_ref_EAF: bool = False,
        temp_dir: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Conditional analysis. Update the beta, se, pval of the conditional SNPs.

        Using cojo-cond in GCTA.
        Condition on the SNPs in cond_snps.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        ldref : str
            LD reference file.
        cond_snps : List[str]
            SNPs to condition on.
        out : str
            Output file.
        sample_size : int
            Sample size.
        use_ref_EAF : bool, optional
            Use the EAF in the LD reference file, by default False

        Returns
        -------
        pd.DataFrame
            Conditional summary statistics.
        """
        if not use_ref_EAF and ColName.EAF not in sumstats.columns:
            raise ValueError(f"{ColName.EAF} is not in the sumstats, please set use_ref_EAF to True")
        chrom = sumstats[ColName.CHR].iloc[0]
        ld = LDRef()
        cojo_input = ld.intersect(sumstats, ldref, f"{temp_dir}/cojo_input_{chrom}", use_ref_EAF)
        cojo_input[ColName.N] = sample_size
        cojo_input = cojo_input[
            [ColName.SNPID, ColName.EA, ColName.NEA, ColName.EAF, ColName.BETA, ColName.SE, ColName.P, ColName.N]
        ]
        cojo_input.rename(
            columns={
                ColName.SNPID: "SNP",
                ColName.EA: "A1",
                ColName.NEA: "A2",
                ColName.EAF: "freq",
                ColName.BETA: "b",
                ColName.SE: "se",
                ColName.P: "p",
                ColName.N: "N",
            },
            inplace=True,
        )
        cojo_p_file = f"{temp_dir}/cojo_input_{chrom}.ma"
        cojo_input.to_csv(cojo_p_file, sep=" ", index=False)
        with open(f"{temp_dir}/cojo_cond_{chrom}.snps", "w") as f:
            f.write('\n'.join(cond_snps))
        cojo_outfile = f"{temp_dir}/cojo_{chrom}.cond"
        cmd = [
            self.gcta,
            "--bfile",
            ldref,
            "--cojo-file",
            cojo_p_file,
            "--cojo-cond",
            f"{temp_dir}/cojo_cond_{chrom}.snps",
            "--out",
            cojo_outfile,
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            cond_res = pd.read_csv(f"{cojo_outfile}.cond.cma.cojo", sep="\t")
            return cond_res

    @io_in_tempdir('./tmp/finemapping')
    def make_ld(self, sumstat: pd.DataFrame, ldref: LDRef, out: str):
        """Make the LD matrix."""
        self.logger.info("Making the LD matrix")

    def run_finemap(self, sumstat: pd.DataFrame, ldref: LDRef, out: str):
        """Run FINEMAP."""
        self.logger.info("Running FINEMAP")
        raise NotImplementedError

    def run_paintor(self, sumstat: pd.DataFrame, ldref: LDRef, out: str):
        """Run PAINTOR."""
        self.logger.info("Running PAINTOR")
        raise NotImplementedError

    def finemap_locus(
        self,
        locus: pd.DataFrame,
        sumstat: pd.DataFrame,
        ldref: str,
        out: str,
        method: str,
        sample_size: int,
        use_ref_EAF: bool = False,
        temp_dir: Optional[str] = None,
    ) -> None:
        """
        Finemap a locus.

        Parameters
        ----------
        locus : Loci
            Locus to finemap.
        sumstat : SumStat
            Summary statistics.
        ldref : LDRef
            LD reference.
        out : str
            Output file.
        method : str
            Method to use.
        sample_size : int
            Sample size.
        use_ref_EAF : bool, optional
            Use the EAF in the LD reference file, by default False
        temp_dir : Optional[str], optional
            Temporary directory, by default None

        Returns
        -------
        pd.DataFrame
            Finemapping results.
        """
        if method not in ["finemap", "paintor", "caviarbf", "susie", "polyfun+susie"]:
            raise ValueError(f"Method {method} not supported")
        if method == "finemap":
            raise NotImplementedError
        elif method == "paintor":
            raise NotImplementedError
        elif method == "caviarbf":
            raise NotImplementedError
        elif method == "susie":
            raise NotImplementedError
        elif method == "polyfun+susie":
            raise NotImplementedError
