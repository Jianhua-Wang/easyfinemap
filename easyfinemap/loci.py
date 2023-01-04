"""Get the independent loci from the input file.

support three approaches:
1. identify the independent lead snps by distance only,
TODO:2. identify the independent lead snps by LD clumping,
TODO:3. identify the independent lead snps by conditional analysis.

than expand the independent lead snps to independent loci by given range.
merge the overlapped independent loci (optional).
"""

import logging
import tempfile
from pathlib import Path
from subprocess import PIPE, run
from typing import List, Optional, Union

import pandas as pd

from easyfinemap.constant import ColName
from easyfinemap.tools import Tools
from easyfinemap.utils import io_in_tempdir, make_SNPID_unique


class Loci:
    """Identify the independent loci."""

    def __init__(self):
        """Initialize the Loci class."""
        self.logger = logging.getLogger("Loci")
        self.plink = Tools().plink
        self.tmp_root = Path.cwd() / "tmp" / "loci"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)

    def identify_indep_loci(
        self,
        sig_df: pd.DataFrame,
        method: str = "distance",
        distance: int = 1000000,
        range: int = 1000000,
        merge: bool = True,
    ) -> pd.DataFrame:
        """
        Identify the independent loci.

        Parameters
        ----------
        sig_df : pd.DataFrame
            The significant snps.
        method : str, optional
            The method to identify the independent loci, by default "distance"
        """
        raise NotImplementedError

    @staticmethod
    def merge_overlapped_loci(loci_df: pd.DataFrame):
        """
        Merge the overlapped loci.

        Parameters
        ----------
        loci_df : pd.DataFrame
            The independent loci.

        Returns
        -------
        pd.DataFrame
            The merged independent loci.
        """
        merged_loci = loci_df.copy()
        merged_loci.sort_values([ColName.CHR, ColName.START, ColName.END], inplace=True)
        merged_loci['no_overlap'] = merged_loci[ColName.START] > merged_loci[ColName.END].shift().cummax()
        merged_loci['diff_chr'] = merged_loci[ColName.CHR] != merged_loci[ColName.CHR].shift()
        merged_loci["break"] = merged_loci["no_overlap"] | merged_loci['diff_chr']
        merged_loci['group'] = merged_loci['break'].cumsum()
        merged_loci = merged_loci.sort_values(['group', ColName.LEAD_SNP_P], ascending=True)
        agg_func = {}
        for col in loci_df.columns:
            if col == ColName.START:
                agg_func[col] = 'min'
            elif col == ColName.END:
                agg_func[col] = 'max'
            else:
                agg_func[col] = 'first'
        result = merged_loci.groupby("group").agg(agg_func)
        result.reset_index(drop=True, inplace=True)
        return result

    @staticmethod
    def indep_snps_by_distance(sig_df: pd.DataFrame, distance: int = 500000) -> pd.DataFrame:
        """
        Identify the independent snps by distance only.

        Parameters
        ----------
        sig_df : pd.DataFrame
            The significant snps.
        distance : int, optional
            The distance threshold, by default 1000000

        Returns
        -------
        pd.DataFrame
            The independent snps.
        """
        sig_df.sort_values(ColName.P, inplace=True)
        lead_snp = []
        while len(sig_df):
            lead_snp.append(sig_df.iloc[[0]])
            sig_df = sig_df[
                ~(
                    (sig_df[ColName.CHR] == sig_df.iloc[0][ColName.CHR])
                    & (sig_df[ColName.BP] >= sig_df.iloc[0][ColName.BP] - distance)
                    & (sig_df[ColName.BP] <= sig_df.iloc[0][ColName.BP] + distance)
                )
            ]  # type: ignore
        lead_snp = pd.concat(lead_snp, axis=0, ignore_index=True)
        return lead_snp

    @staticmethod
    def indep_snps_by_ldclumping(
        sig_df: pd.DataFrame, ldref: str, clump_p1: float = 5e-8, clump_kb: int = 500000, clump_r2: float = 0.1
    ) -> pd.DataFrame:
        """
        Identify the independent snps by LD clumping.

        Parameters
        ----------
        sig_df : pd.DataFrame
            The significant snps.
        ldref : str
            The LD reference file, (plink bfile format, containing wildcard {chrom}), e.g. EUR.chr{chrom}.
        clump_p1 : float, optional
            The p1 threshold, by default 5e-8
        clump_kb : int, optional
            The kb threshold, by default 500000
        clump_r2 : float, optional
            The r2 threshold, by default 0.1

        Returns
        -------
        pd.DataFrame
        """
        clumped_snps = []
        for chrom in sig_df[ColName.CHR].unique():
            sig_df_chr = sig_df[sig_df[ColName.CHR] == chrom]
            clumped_snps.append(Loci().clump_per_chr(sig_df_chr, ldref, clump_p1, clump_kb, clump_r2))  # type: ignore
        clumped_snps = pd.concat(clumped_snps, axis=0, ignore_index=True)
        return clumped_snps

    @io_in_tempdir(dir="./tmp/loci")
    def clump_per_chr(
        self,
        sig_df: pd.DataFrame,
        ldref: str,
        clump_p1: float,
        clump_kb: int,
        clump_r2: float,
        temp_dir: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        LD clumping per chromosome.

        Parameters
        ----------
        sig_df : pd.DataFrame
            The significant snps.
        ldref : str
            The LD reference file, (plink bfile format, containing wildcard {chrom}), e.g. EUR.chr{chrom}.
        clump_p1 : float
            The p1 threshold.
        clump_kb : int
            The kb threshold.
        clump_r2 : float
            The r2 threshold.
        temp_dir : Optional[str], optional
            The temporary directory, by default None

        Returns
        -------
        pd.DataFrame
            The clumped snps.
        """
        chrom = sig_df[ColName.CHR].unique()[0]
        clump_p_file = f"{temp_dir}/clump_p_{chrom}.txt"
        sig_df[[ColName.SNPID, ColName.P]].to_csv(clump_p_file, sep="\t", index=False)
        clump_outfile = f"{temp_dir}/clump_{chrom}.clumped"
        cmd = [
            self.plink,
            "--bfile",
            ldref.format(chrom=chrom),
            "--clump",
            clump_p_file,
            "--clump-p1",
            str(clump_p1),
            "--clump-kb",
            str(clump_kb),
            "--clump-r2",
            str(clump_r2),
            "--clump-snp-field",
            ColName.SNPID,
            "--clump-field",
            ColName.P,
            "--out",
            f"{temp_dir}/clump_{chrom}",
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            clump_snps = pd.read_csv(clump_outfile, delim_whitespace=True, usecols=["SNP"])
            clump_snps = clump_snps["SNP"].to_list()
            clump_snps = sig_df[sig_df[ColName.SNPID].isin(clump_snps)]
            return clump_snps

    @staticmethod
    def indep_snps_by_conditional(sig_df: pd.DataFrame, ld_df: pd.DataFrame, r2_threshold: float = 0.8) -> pd.DataFrame:
        """Identify the independent snps by conditional analysis."""
        raise NotImplementedError

    def cojo_per_chr(self, sig_df: pd.DataFrame, ldref: str, temp_dir: Optional[str] = None) -> pd.DataFrame:
        """
        Conditional analysis per chromosome.

        Parameters
        ----------
        sig_df : pd.DataFrame
            The significant snps.
        ldref : str
            The LD reference file, (plink bfile format, containing wildcard {chrom}), e.g. EUR.chr{chrom}.
        temp_dir : Optional[str], optional
            The temporary directory, by default None

        Returns
        -------
        pd.DataFrame
            The conditional snps.
        """
        chrom = sig_df[ColName.CHR].unique()[0]
        cojo_p_file = f"{temp_dir}/cojo_p_{chrom}.txt"
        sig_df[[ColName.SNPID, ColName.P]].to_csv(cojo_p_file, sep="\t", index=False)
        cojo_outfile = f"{temp_dir}/cojo_{chrom}.cojo"
        cmd = [
            self.plink,
            "--bfile",
            ldref.format(chrom=chrom),
            "--cojo-file",
            cojo_p_file,
            "--cojo-slct",
            "--cojo-p",
            "1",
            "--out",
            f"{temp_dir}/cojo_{chrom}",
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            cojo_snps = pd.read_csv(cojo_outfile, delim_whitespace=True, usecols=["SNP"])
            cojo_snps = cojo_snps["SNP"].to_list()
            cojo_snps = sig_df[sig_df[ColName.SNPID].isin(cojo_snps)]
            return cojo_snps

    @staticmethod
    def leadsnp2loci(sig_df: pd.DataFrame, range: int = 500000, if_merge: bool = True) -> pd.DataFrame:
        """
        Expand the independent lead snps to independent loci by given range.

        Parameters
        ----------
        sig_df : pd.DataFrame
            The independent lead snps.
        range : int, optional
            The range, by default 500000
        if_merge : bool, optional
            Whether merge the overlapped loci, by default True

        Returns
        -------
        pd.DataFrame
            The independent loci.
        """
        loci_df = sig_df.copy()
        loci_df = make_SNPID_unique(loci_df)
        loci_df = loci_df[[ColName.CHR, ColName.BP, ColName.P, ColName.SNPID]]
        loci_df.columns = [ColName.CHR, ColName.LEAD_SNP_BP, ColName.LEAD_SNP_P, ColName.LEAD_SNP]  # type: ignore
        loci_df[ColName.START] = loci_df[ColName.LEAD_SNP_BP] - range
        loci_df[ColName.START] = loci_df[ColName.START].apply(lambda x: 0 if x < 0 else x)
        loci_df[ColName.END] = loci_df[ColName.LEAD_SNP_BP] + range
        loci_df = loci_df[ColName.loci_cols].copy()
        if if_merge:
            loci_df = Loci.merge_overlapped_loci(loci_df)
        return loci_df
