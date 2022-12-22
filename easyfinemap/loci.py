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
from easyfinemap.constant import ColName

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
    sig_df["locus"] = (sig_df["chr"] * 1e9 + sig_df["pos"] // distance).astype(int)
    return sig_df

def indep_snps_by_ld(sig_df: pd.DataFrame, ld_df: pd.DataFrame, r2_threshold: float = 0.8) -> pd.DataFrame:
    raise NotImplementedError

def expand_loci(sig_df: pd.DataFrame, range: int = 1000000) -> pd.DataFrame:
    """
    Expand the independent lead snps to independent loci by given range.

    Parameters
    ----------
    sig_df : pd.DataFrame
        The independent lead snps.
    range : int, optional
        The range, by default 1000000

    Returns
    -------
    pd.DataFrame
        The independent loci.
    """
    loci_df = sig_df.copy()
    loci_df[ColName.START] = loci_df[ColName.LEAD_SNP_BP] - range
    loci_df[ColName.START] = loci_df[ColName.START].apply(lambda x: 0 if x < 0 else x)
    loci_df[ColName.END] = loci_df[ColName.LEAD_SNP_BP] + range
    return loci_df
