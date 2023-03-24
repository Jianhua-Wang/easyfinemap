"""Prepare LD reference for easyfinemap.

1. validate the LD reference.
    1.1. remove duplicate SNPs.
    1.2. make SNP names unique, chr-bp-sorted(EA,NEA).
2. intersect the significant snps with the LD reference.
3. make a plink file from the intersected snps.
TODO: 4. calculate LD matrix from the plink file.
5. perform conditional analysis in conditional mode.
"""

import logging
import os
import shutil
from pathlib import Path
from subprocess import PIPE, check_output, run
from typing import List, Optional, Union

import numpy as np
import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool

from easyfinemap.constant import CHROMS, ColName
from easyfinemap.tools import Tools
from easyfinemap.utils import io_in_tempdir, make_SNPID_unique


class LDRef:
    """Prepare LD reference for easyfinemap."""

    def __init__(self):
        """Initialize the LDRef class."""
        self.logger = logging.getLogger("LDRef")
        self.plink = Tools().plink
        self.gcta = Tools().gcta
        self.tmp_root = Path.cwd() / "tmp" / "ldref"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)

    @io_in_tempdir(dir='./tmp/ldref')
    def _clean_per_chr(self, inprefix: str, outprefix: str, mac: int = 10, temp_dir: Optional[str] = None) -> None:
        """
        Clean the extracted LD reference per chromosome.

        1. Remove duplicated snps.
        2. Make SNP names unique, chr-bp-sorted(EA,NEA).

        Parameters
        ----------
        prefix : str
            The prefix of the extracted LD reference.
        outprefix : str, optional
            The prefix of the cleaned LD reference.
        mac : int, optional
            The minor allele count threshold, by default 10
            SNPs with MAC < mac will be removed.
        temp_dir : Optional[str], optional
            The temp dir, by default None

        Returns
        -------
        None
        """
        prefix = f"{temp_dir}/{inprefix.split('/')[-1]}"
        bim_file = f"{inprefix}.bim"
        bim = pd.read_csv(
            bim_file,
            delim_whitespace=True,
            names=[ColName.CHR, ColName.RSID, "cM", ColName.BP, ColName.EA, ColName.NEA],
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
        shutil.move(f"{prefix}.clean.bed", f"{outprefix}.bed")
        shutil.move(f"{prefix}.clean.bim", f"{outprefix}.bim")
        shutil.move(f"{prefix}.clean.fam", f"{outprefix}.fam")

    @io_in_tempdir(dir='./tmp/ldref')
    def valid(
        self,
        ldref_path: str,
        outprefix: str,
        file_type: str = "plink",
        mac: int = 10,
        threads: int = 1,
        temp_dir: Optional[str] = None,
    ) -> None:
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
            SNPs with MAC < mac will be removed.
        threads : int, optional
            The number of threads to use, by default 1
        temp_dir : Optional[str], optional
            The path to the temporary directory, by default None

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
                    if chrom == 23:
                        inprefix = ldref_path.replace("{chrom}", "X")
                        if os.path.exists(f"{inprefix}.bed"):
                            self.logger.warning(f"chr{chrom} not found, use X instead.")
                            params[0].append(inprefix)
                            params[1].append(f"{outprefix}.chr{chrom}")
                            params[2].append(mac)
                        else:
                            self.logger.warning(f"{inprefix}.bed not found.")
                    else:
                        self.logger.warning(f"{inprefix}.bed not found.")
                        continue
                else:
                    params[0].append(inprefix)
                    params[1].append(f"{outprefix}.chr{chrom}")
                    params[2].append(mac)
            else:
                inprefix = ldref_path
                if not os.path.exists(f"{inprefix}.bed"):
                    raise FileNotFoundError(f"{inprefix}.bed not found.")
                else:
                    # check if chrom is in the bim file
                    res = check_output(f'grep "^{chrom}[[:space:]]" {inprefix}.bim | head -n 1', shell=True)
                    if len(res.decode()) == 0:
                        self.logger.warning(f"Chrom {chrom} not found in {inprefix}.bim")
                        continue
                    else:
                        intermed_prefix = f"{temp_dir}/{outprefix.split('/')[-1]}.chr{chrom}"
                        self.extract(inprefix, intermed_prefix, chrom, mac=mac)
                        params[0].append(intermed_prefix)
                        params[1].append(f"{outprefix}.chr{chrom}")
                        params[2].append(mac)

        with Pool(threads) as p:
            p.map(self._clean_per_chr, *params)

    @io_in_tempdir(dir="./tmp/ldref")
    def extract(
        self,
        inprefix: str,
        outprefix: str,
        chrom: int,
        temp_dir: Optional[str] = None,
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
        temp_dir : str
            The temporary directory.
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
        region_file = f"{temp_dir}/{outprefix.split('/')[-1]}.region"
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

    @io_in_tempdir(dir="./tmp/ldref")
    def intersect(
        self,
        sumstats: pd.DataFrame,
        ldref: str,
        out_plink: str,
        use_ref_EAF: bool = False,
        temp_dir: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Intersect the significant snps with the LD reference.

        Parameters
        ----------
        sumstats : pd.DataFrame
            The summary statistics.
        ldref : str
            The path to the LD reference file.
        out_plink : str
            The output prefix.
        use_ref_EAF : bool, optional
            Use the EAF in the LD reference, by default False
        temp_dir : Optional[str], optional
            The path to the temporary directory, by default None

        Returns
        -------
        pd.DataFrame
            The intersected significant snps.
        """
        if not os.path.exists(f"{ldref}.bim"):
            raise FileNotFoundError(f"{ldref}.bim not found.")
        sumstats[ColName.SNPID].to_csv(f"{temp_dir}/overlap_snpid.txt", index=False, header=False)
        cmd = [
            self.plink,
            "--bfile",
            ldref,
            "--extract",
            f"{temp_dir}/overlap_snpid.txt",
            "--keep-allele-order",
            "--make-bed",
            "--out",
            out_plink,
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(' '.join(cmd))
        self.logger.debug(f"intersect {sumstats.shape[0]} SNPs with {ldref}")
        if res.returncode != 0:
            self.logger.error(res.stderr)
            self.logger.error(f'see log file: {out_plink}.log for details')
            raise RuntimeError(res.stderr)
        bim = pd.read_csv(
            f"{out_plink}.bim",
            delim_whitespace=True,
            names=[ColName.CHR, ColName.RSID, "cM", ColName.BP, ColName.EA, ColName.NEA],
        )
        overlap_sumstat = sumstats[sumstats[ColName.SNPID].isin(bim[ColName.RSID])].copy()
        overlap_sumstat.reset_index(drop=True, inplace=True)

        if use_ref_EAF:
            cmd = [
                self.plink,
                "--bfile",
                out_plink,
                "--freq",
                "--out",
                f"{temp_dir}/freq",
            ]
            res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
            self.logger.debug(f"calculate EAF of {out_plink}")
            self.logger.debug(f"calculate EAF: {' '.join(cmd)}")
            # if res.returncode != 0:
            #     self.logger.error(res.stderr)
            #     self.logger.error(f'see log file: {temp_dir}/freq.log for details')
            #     raise RuntimeError(res.stderr)
            freq = pd.read_csv(f"{temp_dir}/freq.frq", delim_whitespace=True)
            freq['A2_frq'] = 1 - freq['MAF']
            overlap_sumstat['EAF'] = freq['A2_frq'].where(freq['A2'] == overlap_sumstat['EA'], freq['MAF'])
            overlap_sumstat['MAF'] = freq['MAF']
        return overlap_sumstat

    @io_in_tempdir('./tmp/ldref')
    def make_ld(
        self,
        ldref: str,
        outprefix: str,
        **kwargs,
    ):
        """
        Make the LD matrix.

        TODO: Calculate LD matrix using plink-pandas, because plink1.9 --ld contains bug.

        Parameters
        ----------
        ldref : str
            The path to the LD reference file.
        outprefix : str
            The output prefix.

        Raises
        ------
        RuntimeError
            If the return code is not 0.

        Returns
        -------
        None
        """
        self.logger.info("Making the LD matrix")
        cmd = [
            self.plink,
            "--bfile",
            ldref,
            "--r2",
            "square",
            "spaces",
            "--threads",
            "1",
            "--out",
            outprefix,
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(f"get LD matrix: {' '.join(cmd)}")
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            self.logger.debug("LD matrix is made")
            run(["sed", "-i", "s/nan/1e-6/g", f"{outprefix}.ld"])
            # matrix = np.loadtxt(f"{outprefix}.ld")
            # matrix[np.isnan(matrix)] = 1e-6
            # np.savetxt(f"{outprefix}.ld", matrix)

    @io_in_tempdir('./tmp/ldref')
    def cojo_cond(
        self,
        sumstats: pd.DataFrame,
        cond_snps: pd.DataFrame,
        ldref: str,
        sample_size: int,
        use_ref_EAF: bool = False,
        temp_dir: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Conditional analysis. Update the beta, se, pval of the conditional SNPs.

        Parameters
        ----------
        sumstats : pd.DataFrame
            The summary statistics.
        cond_snps : pd.DataFrame
            The conditional SNPs.
        ldref : str
            The path to the LD reference file.
        sample_size : int
            The sample size.
        use_ref_EAF : bool, optional
            Whether to use the EAF in the LD reference file, by default False
        temp_dir : Optional[str], optional
            The path to the temporary directory, by default None

        Raises
        ------
        ValueError
            If the EAF is not in the sumstats and use_ref_EAF is False.

        Returns
        -------
        pd.DataFrame
            The updated summary statistics.
        """
        if not use_ref_EAF and ColName.EAF not in sumstats.columns:
            raise ValueError(f"{ColName.EAF} is not in the sumstats, please set use_ref_EAF to True")
        chrom = sumstats[ColName.CHR].iloc[0]
        ld = LDRef()
        all_sumstats = pd.concat([sumstats, cond_snps], ignore_index=True)
        all_sumstats.drop_duplicates(subset=[ColName.SNPID], inplace=True)
        all_sumstats.sort_values(by=[ColName.CHR, ColName.BP], inplace=True)
        all_sumstats.reset_index(drop=True, inplace=True)
        cojo_input = ld.intersect(all_sumstats, ldref, f"{temp_dir}/cojo_input_{chrom}", use_ref_EAF)
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
            f.write('\n'.join(cond_snps[ColName.SNPID].tolist()))
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
            cond_res = pd.read_csv(f"{cojo_outfile}.cma.cojo", sep="\t", usecols=["SNP", "bC", "bC_se", "pC"])
            cond_res.rename(
                columns={"SNP": ColName.SNPID, "bC": ColName.COJO_BETA, "bC_se": ColName.COJO_SE, "pC": ColName.COJO_P},
                inplace=True,
            )
            output = sumstats.merge(cond_res, on=ColName.SNPID, how="left")
            output = output.dropna(subset=[ColName.COJO_P, ColName.COJO_BETA, ColName.COJO_SE])
            return output
