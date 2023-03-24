"""Get the independent loci from the input file.

support three approaches:
1. identify the independent lead snps by distance only,
2. identify the independent lead snps by LD clumping,
3. identify the independent lead snps by conditional analysis.

than expand the independent lead snps to independent loci by given range.
merge the overlapped independent loci (optional).
"""

import logging
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from subprocess import PIPE, run
from typing import List, Optional, Tuple, Union

import pandas as pd
from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeElapsedColumn

from easyfinemap.constant import ColName
from easyfinemap.ldref import LDRef

# from easyfinemap.sumstat import SumStat
from easyfinemap.tools import Tools
from easyfinemap.utils import get_significant_snps, io_in_tempdir, make_SNPID_unique


class Loci:
    """Identify the independent loci."""

    def __init__(self):
        """Initialize the Loci class."""
        self.logger = logging.getLogger("Loci")
        self.plink = Tools().plink
        self.gcta = Tools().gcta
        self.tmp_root = Path.cwd() / "tmp" / "loci"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)

    def identify_indep_loci(
        self,
        sumstats: pd.DataFrame,
        sig_threshold: float = 5e-8,
        loci_extend: int = 500,
        ldblock: Optional[str] = None,
        if_merge: bool = False,
        outprefix: Optional[str] = None,
        ldref: Optional[str] = None,
        method: str = "distance",
        distance: int = 500,
        clump_kb: int = 500,
        clump_r2: float = 0.1,
        sample_size: Optional[int] = None,
        cojo_window_kb: int = 10000,
        cojo_collinear: float = 0.9,
        diff_freq: float = 0.2,
        only_use_sig_snps: bool = False,
        use_ref_EAF: bool = False,
        threads: int = 1,
    ) -> Union[Tuple[pd.DataFrame, pd.DataFrame], None]:
        """
        Identify the independent loci.

        Parameters
        ----------
        sumstats : pd.DataFrame
            The input summary statistics.
        sig_threshold : float, optional
            The pvalue threshold, by default 5e-8
        loci_extend : int, optional
            The range to extend the independent lead snps to independent loci, by default 500, unit: kb
        if_merge : bool, optional
            Whether to merge the overlapped independent loci, by default False
        ldref : Optional[str], optional
            The LD reference file, by default None
        method : str, optional
            The method to identify the independent loci, by default "distance",
            choose from ["distance", "clumping", "conditional"]
        distance : int, optional
            The distance threshold to identify the independent loci, by default 500, unit: kb
        clump_kb : int, optional
            The distance threshold for LD clumping, by default 10000, unit: kb
        clump_r2 : float, optional
            The r2 threshold for LD clumping, by default 0.1
        sample_size : Optional[int], optional
            The sample size for conditional analysis, by default None
        cojo_window_kb : int, optional
            The distance threshold for conditional analysis, by default 10000
        cojo_collinear : float, optional
            The collinear threshold for conditional analysis, by default 0.9
        diff_freq : float, optional
            The difference frequency threshold for conditional analysis, by default 0.2
        only_use_sig_snps : bool, optional
            Whether to use the significant snps for conditional analysis, by default False
        use_ref_EAF : bool, optional
            Whether to use the reference EAF for conditional analysis, by default False
        threads : int, optional
            The number of threads, by default 1

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            The independent lead snps and independent loci.
        """
        sumstats = make_SNPID_unique(sumstats)
        if ldblock is not None:
            ldblock = pd.read_csv(ldblock, sep="\t", names=[ColName.CHR, ColName.START, ColName.END])
        if method == "distance":
            sig_df = get_significant_snps(sumstats, sig_threshold)
            lead_snp = self.indep_snps_by_distance(sig_df, distance, ldblock)
        elif method == "clumping":
            clump_p1 = sig_threshold
            if ldref is not None:
                sig_df = get_significant_snps(sumstats, sig_threshold)
                lead_snp = self.indep_snps_by_ldclumping(sig_df, ldref, clump_p1, clump_kb, clump_r2)
            else:
                raise ValueError(f"Please provide the ldref file for method: {method}")
        elif method == "conditional":
            if ldref is None:
                raise ValueError("Please provide the ldref file for conditional analysis.")
            if sample_size is None:
                raise ValueError("Please provide the sample size for conditional analysis.")
            else:
                lead_snp = self.indep_snps_by_conditional(
                    sumstats,
                    ldref,
                    sample_size,
                    sig_threshold,
                    cojo_window_kb,
                    cojo_collinear,
                    diff_freq,
                    use_ref_EAF,
                    only_use_sig_snps,
                    ldblock,
                    threads,
                )
        else:
            raise ValueError(f"Unsupported method: {method}")
        loci = self.leadsnp2loci(lead_snp, loci_extend, if_merge, ldblock)
        if if_merge and ColName.COJO_BETA in lead_snp.columns:
            self.logger.warning("The loci identified by cojo may not need merge.")
            lead_snp = lead_snp[lead_snp[ColName.SNPID].isin(loci[ColName.LEAD_SNP])]
        if outprefix:
            loci_file = f"{outprefix}.loci.txt"
            loci.to_csv(loci_file, sep="\t", index=False)
            self.logger.info(f"Save the independent loci to {loci_file}")
            leadsnp_file = f"{outprefix}.leadsnp.txt"
            lead_snp.to_csv(leadsnp_file, sep="\t", index=False)
            self.logger.info(f"Save the independent lead snps to {leadsnp_file}")
        return lead_snp, loci

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
    def indep_snps_by_distance(
        sig_df: pd.DataFrame, distance: int = 500, ldblock: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Identify the independent snps by distance only.

        Parameters
        ----------
        sig_df : pd.DataFrame
            The significant snps.
        distance : int, optional
            The distance threshold, by default 500, unit: kb
        ldblock : Optional[pd.DataFrame], optional
            The ld block information, use boundary to identify the independent snps, by default None

        Returns
        -------
        pd.DataFrame
            The independent snps.
        """
        sig_df = sig_df.sort_values(ColName.P).copy()
        lead_snp = []
        if ldblock is not None:
            while len(sig_df):
                lead_snp.append(sig_df.iloc[[0]])
                sig_block = ldblock[
                    (ldblock[ColName.CHR] == sig_df.iloc[0][ColName.CHR])
                    & (ldblock[ColName.START] <= sig_df.iloc[0][ColName.BP])
                    & (ldblock[ColName.END] >= sig_df.iloc[0][ColName.BP])
                ]
                sig_df = sig_df[
                    ~(
                        (sig_df[ColName.CHR] == sig_df.iloc[0][ColName.CHR])
                        & (sig_df[ColName.BP] >= sig_df.iloc[0][ColName.BP] - sig_block.iloc[0][ColName.START])
                        & (sig_df[ColName.BP] <= sig_df.iloc[0][ColName.BP] + sig_block.iloc[0][ColName.END])
                    )
                ]
        else:
            distance = distance * 1000
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
        sig_df: pd.DataFrame, ldref: str, clump_p1: float = 5e-8, clump_kb: int = 500, clump_r2: float = 0.1
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
            The kb threshold, by default 500, unit: kb
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
            if os.path.exists(clump_outfile):
                clump_snps = pd.read_csv(clump_outfile, delim_whitespace=True, usecols=["SNP"])
                clump_snps = clump_snps["SNP"].to_list()
                clump_snps = sig_df[sig_df[ColName.SNPID].isin(clump_snps)]
                return clump_snps
            else:
                logging.warning(f"No clumped snps found for chromosome {chrom}")
                return pd.DataFrame()

    @staticmethod
    def indep_snps_by_conditional(
        sumstats: pd.DataFrame,
        ldref: str,
        sample_size: int,
        sig_threshold: float = 5e-8,
        cojo_window_kb: int = 10000,
        cojo_collinear: float = 0.9,
        diff_freq: float = 0.2,
        use_ref_EAF: bool = False,
        only_use_sig_snps: bool = False,
        ldblock: Optional[pd.DataFrame] = None,
        threads: int = 1,
    ) -> pd.DataFrame:
        """
        Identify the independent snps by conditional analysis.

        Parameters
        ----------
        sumstats : pd.DataFrame
            The summary statistics.
        ldref : str
            The LD reference file, (plink bfile format, containing wildcard {chrom}), e.g. EUR.chr{chrom}.
        sample_size : int
            The sample size.
        sig_threshold : float, optional
            The significance threshold, by default 5e-8
        cojo_window_kb : int, optional
            The cojo window, by default 10000, in kb
        cojo_collinear : float, optional
            The cojo collinear, by default 0.9
        diff_freq : float, optional
            The difference frequency, by default 0.2
        use_ref_EAF : bool, optional
            Whether to use the reference EAF, by default False
        only_use_sig_snps : bool, optional
            Whether to only use the significant snps, by default False
        ldblock : Optional[pd.DataFrame], optional
            The LD block, run cojo in each LD block, by default None
        threads : int, optional
            The number of threads, by default 1
        """
        logger = logging.getLogger('COJO')
        if not use_ref_EAF and ColName.EAF not in sumstats.columns:
            raise ValueError(f"{ColName.EAF} is not in the sumstats, please set use_ref_EAF to True")
        sig_df = sumstats[sumstats[ColName.P] <= sig_threshold]
        logger.debug(f"Number of significant snps: {len(sig_df)}")
        logger.debug(f"Number of chromosomes: {len(sig_df[ColName.CHR].unique())}")
        args_list = []
        loci = Loci()
        if ldblock is not None:
            sig_blocks = loci.indep_snps_by_distance(sig_df, ldblock=ldblock)
            sig_blocks = loci.leadsnp2loci(sig_blocks, ldblock=ldblock)
            for i in sig_blocks.index:
                if only_use_sig_snps:
                    in_df = sig_df[
                        (sig_df[ColName.CHR] == sig_blocks.loc[i][ColName.CHR])
                        & (sig_df[ColName.BP] >= sig_blocks.loc[i][ColName.START])
                        & (sig_df[ColName.BP] <= sig_blocks.loc[i][ColName.END])
                    ]
                else:
                    in_df = sumstats[
                        (sumstats[ColName.CHR] == sig_blocks.loc[i][ColName.CHR])
                        & (sumstats[ColName.BP] >= sig_blocks.loc[i][ColName.START])
                        & (sumstats[ColName.BP] <= sig_blocks.loc[i][ColName.END])
                    ]
                args_list.append(
                    (
                        in_df,
                        ldref.format(chrom=sig_blocks.loc[i][ColName.CHR]),
                        sample_size,
                        cojo_window_kb,
                        cojo_collinear,
                        diff_freq,
                        sig_threshold,
                        use_ref_EAF,
                    )
                )
        else:
            for chrom in sig_df[ColName.CHR].unique():
                if only_use_sig_snps:
                    in_df = sig_df[sig_df[ColName.CHR] == chrom]
                else:
                    in_df = sumstats[sumstats[ColName.CHR] == chrom]
                args_list.append(
                    (
                        in_df,
                        ldref.format(chrom=chrom),
                        sample_size,
                        cojo_window_kb,
                        cojo_collinear,
                        diff_freq,
                        sig_threshold,
                        use_ref_EAF,
                    )
                )

        with ProcessPoolExecutor(max_workers=threads) as executor:
            results = []
            with Progress(
                TextColumn("{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TimeElapsedColumn(),
                auto_refresh=False,
            ) as progress:
                task = progress.add_task("Run cojo-slct", total=len(args_list))
                for _ in executor.map(loci.cojo_slct, *zip(*args_list)):
                    progress.update(task, advance=1)
                    progress.refresh()
                    results.append(_)
        cojo_snps = pd.concat(results, axis=0, ignore_index=True)
        return cojo_snps

    @io_in_tempdir(dir="./tmp/loci")
    def cojo_slct(
        self,
        sumstats: pd.DataFrame,
        ldref: str,
        sample_size: int,
        cojo_window_kb: int = 10000,
        cojo_collinear: float = 0.9,
        diff_freq: float = 0.2,
        sig_threshold: float = 5e-8,
        use_ref_EAF: bool = False,
        temp_dir: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Conditional analysis for input sumstatistics.

        Parameters
        ----------
        sumstats : pd.DataFrame
            The input sumstatistics, from same chromosome or locus.
        ldref : str
            The LD reference file, (plink bfile format, containing wildcard {chrom}), e.g. EUR.chr{chrom}.
        sample_size : int
            The sample size of the input sumstatistics.
        cojo_window_kb : int, optional
            The cojo window, by default 10000, unit: kb
        cojo_collinear : float, optional
            The cojo collinear, by default 0.9
        diff_freq : float, optional
            The difference frequency, by default 0.2
        sig_threshold : float, optional
            The significance threshold, by default 5e-8
        use_ref_EAF : bool, optional
            Whether to use the reference EAF, by default False
        temp_dir : Optional[str], optional
            The temporary directory, by default None

        Returns
        -------
        pd.DataFrame
            The conditional snps.
        """
        chrom = sumstats[ColName.CHR].unique()[0]
        cojo_input = sumstats.copy()
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
        cojo_outfile = f"{temp_dir}/cojo_{chrom}.slct"
        cmd = [
            self.gcta,
            "--bfile",
            ldref,
            "--cojo-file",
            cojo_p_file,
            "--cojo-slct",
            "--cojo-p",
            str(sig_threshold),
            "--cojo-wind",
            str(cojo_window_kb),
            "--cojo-collinear",
            str(cojo_collinear),
            "--diff-freq",
            str(diff_freq),
            "--out",
            cojo_outfile,
        ]
        self.logger.debug(f"Run cojo-slct: {' '.join(cmd)}")
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            cojo_snps = pd.read_csv(
                f"{cojo_outfile}.jma.cojo", delim_whitespace=True, usecols=["SNP", "bJ", "bJ_se", "pJ"]
            )
            cojo_snps.rename(
                columns={
                    "SNP": ColName.SNPID,
                    "bJ": ColName.COJO_BETA,
                    "bJ_se": ColName.COJO_SE,
                    "pJ": ColName.COJO_P,
                },
                inplace=True,
            )
            cojo_snps = cojo_snps[cojo_snps[ColName.COJO_P] <= sig_threshold]
            cojo_snps = sumstats.merge(cojo_snps, on=ColName.SNPID, how="inner")
            return cojo_snps

    @staticmethod
    def leadsnp2loci(
        lead_snps: pd.DataFrame, range: int = 500, if_merge: bool = False, ldblock: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Expand the independent lead snps to independent loci by given range.

        Parameters
        ----------
        lead_snps : pd.DataFrame
            The independent lead snps.
        range : int, optional
            The range, by default 500, unit: kb
        if_merge : bool, optional
            Whether merge the overlapped loci, by default False
        ldblock : Optional[pd.DataFrame], optional
            The ld block, using LD block to expand the independent loci, by default None

        Returns
        -------
        pd.DataFrame
            The independent loci.
        """
        loci_df = lead_snps.copy()
        loci_df = loci_df[[ColName.CHR, ColName.BP, ColName.P, ColName.SNPID]]
        loci_df.columns = [ColName.CHR, ColName.LEAD_SNP_BP, ColName.LEAD_SNP_P, ColName.LEAD_SNP]  # type: ignore
        if ldblock is not None:
            loci_df[ColName.START] = 0
            loci_df[ColName.END] = 0
            for i in loci_df.index:
                sub_ldblock = ldblock[
                    (ldblock[ColName.CHR] == loci_df.loc[i, ColName.CHR])
                    & (ldblock[ColName.START] <= loci_df.loc[i, ColName.LEAD_SNP_BP])
                    & (ldblock[ColName.END] >= loci_df.loc[i, ColName.LEAD_SNP_BP])
                ]
                if sub_ldblock.empty:
                    continue
                else:
                    loci_df.loc[i, ColName.START] = sub_ldblock.iloc[0][ColName.START]
                    loci_df.loc[i, ColName.END] = sub_ldblock.iloc[0][ColName.END]
        else:
            range = range * 1000

            loci_df[ColName.START] = loci_df[ColName.LEAD_SNP_BP] - range
            loci_df[ColName.START] = loci_df[ColName.START].apply(lambda x: 0 if x < 0 else x)
            loci_df[ColName.END] = loci_df[ColName.LEAD_SNP_BP] + range
        loci_df = loci_df[ColName.loci_cols].copy()
        if if_merge:
            loci_df = Loci.merge_overlapped_loci(loci_df)
        loci_df = loci_df.sort_values(by=[ColName.CHR, ColName.START, ColName.END])
        return loci_df
