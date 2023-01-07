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
from easyfinemap.utils import get_significant_snps, io_in_tempdir


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
    def run_finemap(
        self,
        sumstats: pd.DataFrame,
        ld_matrix: str,
        sample_size: int,
        max_causal: int = 1,
        temp_dir: Optional[str] = None,
    ) -> pd.Series:
        """
        Run FINEMAP.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        ld_matrix : str
            Path to LD matrix.
        sample_size : int
            Sample size.
        max_causal : int, optional
            Maximum number of causal variants, by default 1
        temp_dir : Optional[str], optional
            Path to tempdir, by default None

        Returns
        -------
        pd.Series
            The result of FINEMAP.
        """
        finemap_input = sumstats.copy()
        finemap_input[ColName.MAF] = finemap_input[ColName.MAF].replace(0, 0.00001)
        finemap_input = finemap_input[
            [ColName.SNPID, ColName.CHR, ColName.BP, ColName.EA, ColName.NEA, ColName.MAF, ColName.BETA, ColName.SE]
        ]
        finemap_input.rename(
            columns={
                ColName.SNPID: "rsid",
                ColName.CHR: "chromosome",
                ColName.BP: "position",
                ColName.MAF: "maf",
                ColName.BETA: "beta",
                ColName.SE: "se",
                ColName.EA: "allele1",
                ColName.NEA: "allele2",
            },
            inplace=True,
        )
        finemap_input.to_csv(f"{temp_dir}/finemap.z", sep=" ", index=False)
        with open(f"{temp_dir}/finemap.master", "w") as f:
            master_content = [
                f"{temp_dir}/finemap.z",
                ld_matrix,
                f"{temp_dir}/finemap.snp",
                f"{temp_dir}/finemap.config",
                f"{temp_dir}/finemap.cred",
                f"{temp_dir}/finemap.log",
                str(sample_size),
            ]
            f.write("z;ld;snp;config;cred;log;n_samples\n")
            f.write(";".join(master_content))
        cmd = [
            self.finemap,
            "--sss",
            "--in-files",
            f"{temp_dir}/finemap.master",
            "--n-causal-snps",
            str(max_causal),
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(f"run FINEMAP: {' '.join(cmd)}")
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            finemap_res = pd.read_csv(f"{temp_dir}/finemap.config", sep=" ", usecols=["config", "prob"])
            finemap_res = pd.Series(finemap_res["prob"].values, index=finemap_res["config"].values)  # type: ignore
            return finemap_res

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

        1. Check if LD is needed, abf, susie, polyfun+susie do not need LD.
        2. If LD is needed, intersect the locus with the LD reference and make the LD matrix.
        3. Run the finemapping method.
        4. Get the finemapping results.
        5. Merge the finemapping results with the input sumstats.
        6. Return the credible set or full summary statistics with posterior probabilities.

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
