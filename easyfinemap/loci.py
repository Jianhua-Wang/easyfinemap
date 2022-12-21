"""Get the independent loci from the input file.
support three approaches:
TODO:1. identify the independent lead snps by distance only,
TODO:2. identify the independent lead snps by LD clumping,
TODO:3. identify the independent lead snps by conditional analysis.

TODO:than expand the independent lead snps to independent loci by given range.
merge the overlapped independent loci (optional).
"""

from easyfinemap.logger import logger
from easyfinemap.tools import Tools

import pandas as pd


def get_significant_snps(df: pd.DataFrame, pvalue_col: str, pvalue_threshold: float = 5e-8):
    """
    Get the significant snps from the input file, filter by pvalue.

    Parameters
    ----------
    df : pd.DataFrame
        The input summary statistics.
    pvalue_col : str
        The column name of pvalue.
    pvalue_threshold : float, optional
        The pvalue threshold, by default 5e-8

    Returns
    -------
    pd.DataFrame
        The significant snps, sorted by pvalue.
    """
    sig_df = df.loc[df[pvalue_col] < pvalue_threshold].copy()
    sig_df.sort_values(pvalue_col, inplace=True)
    return sig_df


def merge_overlapped_loci(loci_df: pd.DataFrame, chr_col: str = "chr", start_col: str = "start", end_col: str = "end"):
    """
    Merge the overlapped loci.
    More details: https://stackoverflow.com/questions/57882621/efficient-merge-overlapping-intervals-in-same-pandas-dataframe-with-start-and-fi

    Parameters
    ----------
    loci_df : pd.DataFrame
        The independent loci.
    chr_col : str, optional
        The column name of chromosome, by default "chr"
    start_col : str, optional
        The column name of start position, by default "start"
    end_col : str, optional
        The column name of end position, by default "end"

    Returns
    -------
    pd.DataFrame
        The merged independent loci.
    """
    merged_loci = loci_df.copy()
    merged_loci.sort_values([chr_col, start_col, end_col], inplace=True)
    merged_loci['no_overlap'] = merged_loci[start_col] > merged_loci[end_col].shift().cummax()
    merged_loci['diff_chr'] = merged_loci[chr_col] != merged_loci[chr_col].shift()
    merged_loci["break"] = merged_loci["no_overlap"] | merged_loci['diff_chr']
    merged_loci['group'] = merged_loci['break'].cumsum()
    result = merged_loci.groupby("group").agg({chr_col: 'max', start_col: "min", end_col: "max"})
    result.reset_index(drop=True, inplace=True)
    return result


def indep_snps_by_distance(sig_df: pd.DataFrame, distance: int = 500000):
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
    sig_df["locus"] = (sig_df["chr"] * 1e9 + sig_df["pos"] // distance).astype(int)
    return sig_df

