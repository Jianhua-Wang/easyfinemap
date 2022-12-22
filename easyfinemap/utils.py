"""Utils for easyfinemap."""

import pandas as pd

from easyfinemap.constant import ColName


def get_significant_snps(df: pd.DataFrame, pvalue_threshold: float = 5e-8):
    """
    Get the significant snps from the input file, filter by pvalue.

    Parameters
    ----------
    df : pd.DataFrame
        The input summary statistics.
    pvalue_threshold : float, optional
        The pvalue threshold, by default 5e-8

    Returns
    -------
    pd.DataFrame
        The significant snps, sorted by pvalue.
    """
    sig_df = df.loc[df[ColName.P] < pvalue_threshold].copy()
    sig_df.sort_values(ColName.P, inplace=True)
    sig_df.reset_index(drop=True, inplace=True)
    return sig_df


def make_SNPID_unique(sumstat: pd.DataFrame, replace_rsIDcol: bool = False, remove_duplicates: bool = True):
    """
    Make the SNPID unique.
    The unique SNPID is chr-bp-sorted(EA,NEA)

    Parameters
    ----------
    sumstat : pd.DataFrame
        The input summary statistics.
    replace_rsIDcol : bool, optional
        Whether to replace the rsID column with the unique SNPID, by default False
    remove_duplicates : bool, optional
        Whether to remove the duplicated SNPs, keep the one with smallest P-value, by default True

    Returns
    -------
    pd.DataFrame
        The summary statistics with unique SNPID.
    """
    df = sumstat.copy()
    allele_df = df[[ColName.EA, ColName.NEA]].copy()
    b = allele_df.values
    b.sort(axis=1)
    allele_df[[ColName.EA, ColName.NEA]] = b
    allele_df[ColName.SNPID] = (
        df[ColName.CHR].astype(str)
        + "-"
        + df[ColName.BP].astype(str)
        + "-"
        + allele_df[ColName.EA]
        + "-"
        + allele_df[ColName.NEA]
    )
    if replace_rsIDcol:
        df[ColName.RSID] = allele_df[ColName.SNPID]
    else:
        df.insert(loc=0, column=ColName.SNPID, value=allele_df[ColName.SNPID].values)  # type: ignore
    if remove_duplicates:
        df.sort_values(ColName.P, inplace=True)
        if replace_rsIDcol:
            df.drop_duplicates(subset=[ColName.RSID], keep="first", inplace=True)
        else:
            df.drop_duplicates(subset=[ColName.SNPID], keep="first", inplace=True)
        df.sort_values([ColName.CHR, ColName.BP], inplace=True)
        df.reset_index(drop=True, inplace=True)
    return df


