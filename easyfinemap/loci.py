"""Get the independent loci from the input file.
support three approaches:
1. identify the independent lead snps by distance only,
TODO:2. identify the independent lead snps by LD clumping,
TODO:3. identify the independent lead snps by conditional analysis.

than expand the independent lead snps to independent loci by given range.
merge the overlapped independent loci (optional).
"""

from easyfinemap.logger import logger
from easyfinemap.tools import Tools
from easyfinemap.constant import ColName
from easyfinemap.utils import make_SNPID_unique

import pandas as pd


def merge_overlapped_loci(loci_df: pd.DataFrame):
    """
    Merge the overlapped loci.
    More details: https://stackoverflow.com/questions/57882621/efficient-merge-overlapping-intervals-in-same-pandas-dataframe-with-start-and-fi

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


def indep_snps_by_ld(sig_df: pd.DataFrame, ld_df: pd.DataFrame, r2_threshold: float = 0.8) -> pd.DataFrame:
    raise NotImplementedError


def leadsnp2loci(sig_df: pd.DataFrame, range: int = 500000, if_merge: bool = True) -> pd.DataFrame:
    """
    Expand the independent lead snps to independent loci by given range.

    Parameters
    ----------
    sig_df : pd.DataFrame
        The independent lead snps.
    range : int, optional
        The range, by default 1000000
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
    loci_df.columns = [ColName.CHR, ColName.LEAD_SNP_BP, ColName.LEAD_SNP_P, ColName.LEAD_SNP] # type: ignore
    loci_df[ColName.START] = loci_df[ColName.LEAD_SNP_BP] - range
    loci_df[ColName.START] = loci_df[ColName.START].apply(lambda x: 0 if x < 0 else x)
    loci_df[ColName.END] = loci_df[ColName.LEAD_SNP_BP] + range
    loci_df = loci_df[ColName.loci_cols].copy()
    if if_merge:
        loci_df = merge_overlapped_loci(loci_df)
    return loci_df
