"""Prepare LD reference for easyfinemap.
TODO: 1. intersect the significant snps with the LD reference.
TODO: 2. make a plink file from the intersected snps.
TODO: 3. calculate LD matrix from the plink file.
"""

from subprocess import run, PIPE
from pathlib import Path
import tempfile
import logging

import pandas as pd
import numpy as np

from easyfinemap.constant import ColName
from easyfinemap.tools import Tools
from easyfinemap.utils import get_significant_snps, make_SNPID_unique

import logging


class LDRef:
    """Prepare LD reference for easyfinemap."""

    def __init__(self, ldref_path: str, file_type: str = "plink", log_level: str = "DEBUG"):
        """
        Initialize the LDRef class.

        Parameters
        ----------
        ldref_path : str
            The path to the LD reference file.
        file_type : str, optional
            The file type of the LD reference file, by default "plink"
            TODO: support other file types, e.g. vcf.
        log_level : str, optional
            The log level, by default "DEBUG"
        """
        self.logger = logging.getLogger(f"LDRef")
        self.logger.setLevel(log_level)
        self.ldref_path = ldref_path
        tools = Tools()
        self.plink = tools.plink
        # TODO: support input files with wildcards.
        # if "{chrom}" not in self.ldref_path.name:
        #     raise ValueError("The LD reference file name should contain {chrom}.")
        self.file_type = file_type
        self.tmp_root = Path.cwd() / "tmp" / "ldref"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)
        # self.ldref_dir = Path(tempfile.mkdtemp())
        self.temp_dir = tempfile.mkdtemp(dir=self.tmp_root)
        self.logger.debug(f"LDRef temp dir: {self.temp_dir}")
        self.temp_dir_path = Path(self.temp_dir)

    def clean(self, prefix: str) -> None:
        """
        Clean the extracted LD reference.
        1. Remove duplicated snps.
        2. Make SNP names unique, chr-bp-sorted(EA,NEA).

        Parameters
        ----------
        prefix : str
            The prefix of the extracted LD reference.
        """
        prefix = f"{self.temp_dir_path}/{prefix}"
        bim_file = f"{prefix}.bim"
        bim = pd.read_csv(
            bim_file, sep="\t", names=[ColName.CHR, ColName.RSID, "cM", ColName.BP, ColName.EA, ColName.NEA]
        )
        bim[1] = bim.index  # use number as rsid, make sure it is unique
        bim.to_csv(bim_file, sep="\t", index=False, header=False)
        # TODO: remove duplicated snps.
        cmd = [
            self.plink,
            "--bfile",
            prefix,
            "--keep-allele-order",
            "--list-duplicate-vars",
            "ids-only",
            "suppress-first",
            "--out",
            f"{prefix}",
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(' '.join(cmd))
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        cmd = [
            self.plink,
            "--bfile",
            f"{prefix}",
            "--exclude",
            f"{prefix}.dupvar",
            "--keep-allele-order",
            "--make-bed",
            "--out",
            f"{prefix}.clean",
        ]
        self.logger.debug(' '.join(cmd))
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        rmdup_bim = pd.read_csv(
            f"{prefix}.clean.bim",
            delim_whitespace=True,
            names=[ColName.CHR, ColName.RSID, "cM", ColName.BP, ColName.EA, ColName.NEA],
        )
        rmdup_bim = make_SNPID_unique(rmdup_bim, replace_rsIDcol=True, remove_duplicates=False)
        rmdup_bim.to_csv(f"{prefix}.clean.bim", sep="\t", index=False, header=False)

    def extract(self, chrom: int, start: int, end: int, outprefix: str, mac: int = 10) -> None:
        """
        Extract the genotypes of given region from the LD reference.

        Parameters
        ----------
        chrom : int
            The chromosome number.
        start : int
            The start position.
        end : int
            The end position.
        outprefix : str
            The output prefix.
        mac : int, optional
            The minor allele count threshold, by default 10
        """
        outprefix = f"{self.temp_dir_path}/{outprefix}"
        region_file = f"{outprefix}.region"
        with open(region_file, "w") as f:
            f.write(f"{chrom}\t{start}\t{end}\tregion")
        cmd = [
            self.plink,
            "--bfile",
            self.ldref_path,
            "--extract",
            "range",
            str(region_file),
            "--keep-allele-order",
            "--mac",
            str(mac),
            "--make-bed",
            "--out",
            outprefix,
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(' '.join(cmd))
        self.logger.debug(f"extract {chrom}:{start}-{end} from {self.ldref_path}")
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)

    def intersect(self, sig_snps: pd.DataFrame) -> pd.DataFrame:
        """
        Intersect the significant snps with the LD reference.

        Parameters
        ----------
        sig_snps : pd.DataFrame
            The significant snps.

        Returns
        -------
        pd.DataFrame
            The intersected significant snps.
        """
        raise NotImplementedError
