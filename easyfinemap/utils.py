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
    return sig_df

