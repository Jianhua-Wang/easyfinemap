"""Prepare LD reference for easyfinemap.

1. validate the LD reference.
    1.1. remove duplicate SNPs.
    1.2. make SNP names unique, chr-bp-sorted(EA,NEA).
TODO: 1. intersect the significant snps with the LD reference.
TODO: 2. make a plink file from the intersected snps.
TODO: 3. calculate LD matrix from the plink file.
"""

import logging
import os
import shutil
import tempfile
from pathlib import Path
from subprocess import PIPE, run
from typing import List, Optional, Union

import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool

from easyfinemap.constant import CHROMS, ColName
from easyfinemap.tools import Tools
from easyfinemap.utils import make_SNPID_unique


class LDRef:
    """Prepare LD reference for easyfinemap."""

    def __init__(self, log_level: str = "DEBUG"):
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
        self.logger = logging.getLogger("LDRef")
        self.logger.setLevel(log_level)
        self.plink = Tools().plink
        self.tmp_root = Path.cwd() / "tmp" / "ldref"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)
        self.temp_dir = tempfile.mkdtemp(dir=self.tmp_root)
        self.logger.debug(f"LDRef temp dir: {self.temp_dir}")
        self.temp_dir_path = Path(self.temp_dir)

    def clean(self, inprefix: str, outprefix: Optional[str] = None, mac: int = 10) -> None:
        """
        Clean the extracted LD reference.

        1. Remove duplicated snps.
        2. Make SNP names unique, chr-bp-sorted(EA,NEA).

        Parameters
        ----------
        prefix : str
            The prefix of the extracted LD reference.
        outprefix : str, optional
            The prefix of the cleaned LD reference, by default None
            If None, the cleaned LD reference will be saved to the same directory as the extracted LD reference.
        mac : int, optional
            The minor allele count threshold, by default 10
            SNPs with MAC < mac will be removed.

        Returns
        -------
        None
        """
        prefix = f"{self.temp_dir_path}/{inprefix.split('/')[-1]}"
        bim_file = f"{inprefix}.bim"
        bim = pd.read_csv(
            bim_file, sep="\t", names=[ColName.CHR, ColName.RSID, "cM", ColName.BP, ColName.EA, ColName.NEA]
        )
        bim[ColName.RSID] = bim.index  # use number as rsid, make sure it is unique
        bim.to_csv(f"{prefix}.bim", sep="\t", index=False, header=False)

        if mac <= 0:
            raise ValueError(f"mac should be > 0, got {mac}.")
        cmd = [
            self.plink,
            "--bed",
            f"{inprefix}.bed",
            "--fam",
            f"{inprefix}.fam",
            "--bim",
            f"{prefix}.bim",
            "--keep-allele-order",
            "--list-duplicate-vars",
            "ids-only",
            "suppress-first",
            "--out",
            f"{prefix}",
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(' '.join(cmd))
        cmd = [
            self.plink,
            "--bed",
            f"{inprefix}.bed",
            "--fam",
            f"{inprefix}.fam",
            "--bim",
            f"{prefix}.bim",
            "--exclude",
            f"{prefix}.dupvar",
            "--mac",
            str(mac),
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
        if outprefix is not None:
            shutil.move(f"{prefix}.clean.bed", f"{outprefix}.bed")
            shutil.move(f"{prefix}.clean.bim", f"{outprefix}.bim")
            shutil.move(f"{prefix}.clean.fam", f"{outprefix}.fam")

    def valid(self, ldref_path: str, outprefix: str, file_type: str = "plink", mac: int = 10, threads: int = 1) -> None:
        """
        Validate the LD reference file.

        TODO:1. format vcfs to plink files.
        2. remove duplicated snps.
        3. remove snps with MAC < mac.
        4. make SNP names unique, chr-bp-sorted(EA,NEA).
        TODO:5. mark bim file with "#easyfinemap validated" flag in the first line.

        Parameters
        ----------
        ldref_path : str
            The path to the LD reference file.
        outprefix : str
            The output prefix.
        file_type : str, optional
            The file type of the LD reference file, by default "plink"
        mac: int, optional
            The minor allele count threshold, by default 10

        Raises
        ------
        ValueError
            If the file type is not supported.

        Returns
        -------
        None
        """
        if file_type == "plink":
            self.file_type = file_type
        else:
            raise ValueError(f"Unsupported file type: {file_type}")

        params: List[List[Union[str, int]]] = [[] for _ in range(3)]
        for chrom in CHROMS:
            if "{chrom}" in ldref_path:
                inprefix = ldref_path.replace("{chrom}", str(chrom))
                if not os.path.exists(f"{inprefix}.bed"):
                    raise FileNotFoundError(f"{inprefix}.bed not found.")
                else:
                    params[0].append(inprefix)
                    params[1].append(f"{outprefix}.chr{chrom}")
                    params[2].append(mac)
            else:
                inprefix = ldref_path
                if not os.path.exists(f"{inprefix}.bed"):
                    raise FileNotFoundError(f"{inprefix}.bed not found.")
                else:
                    intermed_prefix = f"{self.temp_dir}/{outprefix.split('/')[-1]}.chr{chrom}"
                    self.extract(inprefix, intermed_prefix, chrom, mac=mac)
                    params[0].append(intermed_prefix)
                    params[1].append(f"{outprefix}.chr{chrom}")
                    params[2].append(mac)

        with Pool(threads) as p:
            p.map(self.clean, *params)

    def extract(
        self,
        inprefix: str,
        outprefix: str,
        chrom: int,
        start: Optional[int] = None,
        end: Optional[int] = None,
        mac: int = 10,
    ) -> None:
        """
        Extract the genotypes of given region from the LD reference.

        Parameters
        ----------
        inprefix : str
            The input prefix.
        outprefix : str
            The output prefix.
        chrom : int
            The chromosome number.
        start : int, optional
            The start position, by default None
        end : int, optional
            The end position, by default None
        mac: int, optional
            The minor allele count threshold, by default 10

        Returns
        -------
        None
        """
        region_file = f"{self.temp_dir}/{outprefix.split('/')[-1]}.region"
        if start is None:
            extract_cmd = ["--chr", str(chrom)]
        else:
            with open(region_file, "w") as f:
                f.write(f"{chrom}\t{start}\t{end}\tregion")
            extract_cmd = ["--extract", "range", region_file]

        if "{chrom}" in inprefix:
            inprefix = inprefix.replace("{chrom}", str(chrom))
        if not os.path.exists(f"{inprefix}.bed"):
            raise FileNotFoundError(f"{inprefix}.bed not found.")
        cmd = [
            self.plink,
            "--bfile",
            inprefix,
            *extract_cmd,
            "--keep-allele-order",
            "--mac",
            str(mac),
            "--make-bed",
            "--out",
            outprefix,
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(' '.join(cmd))
        self.logger.debug(f"extract chr{chrom}:{start}-{end} from {inprefix}")
        if res.returncode != 0:
            self.logger.error(res.stderr)
            self.logger.error(f'see log file: {outprefix}.log for details')
            raise RuntimeError(res.stderr)

    def intersect(self, sig_snps: pd.DataFrame, use_ref_EAF: bool = False) -> pd.DataFrame:
        """
        Intersect the significant snps with the LD reference.

        Parameters
        ----------
        sig_snps : pd.DataFrame
            The significant snps.
        use_ref_EAF : bool, optional
            Use the EAF in the LD reference, by default False

        Returns
        -------
        pd.DataFrame
            The intersected significant snps.
        """
        raise NotImplementedError
