"""Tests for the utils module."""

import pytest

from easyfinemap.utils import get_significant_snps, make_SNPID_unique

def test_get_significant_snps(sumstats_data):
    """Test the get_significant_snps function."""
    sig_df = get_significant_snps(sumstats_data)
    assert len(sig_df) == 9

def test_make_SNPID_unique(mock_sumstat):
    """Test the make_SNPID_unique function."""
    unique_df1 = make_SNPID_unique(mock_sumstat, replace_rsIDcol=False)
    unique_df2 = make_SNPID_unique(mock_sumstat, replace_rsIDcol=True)
    unique_df3 = make_SNPID_unique(mock_sumstat, remove_duplicates=False)
    assert len(unique_df1) == 8
    assert len(unique_df2) == 8
    assert len(unique_df3) == 9
    assert "SNPID" in unique_df1.columns
    assert "SNPID" not in unique_df2.columns