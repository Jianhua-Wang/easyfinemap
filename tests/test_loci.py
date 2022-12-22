"""Tests for the loci module."""

import pytest

from easyfinemap.loci import merge_overlapped_loci, indep_snps_by_distance, leadsnp2loci
from easyfinemap.constant import ColName

def test_merge_overlapped_loci(loci_data):
    """Test the merge_overlapped_loci function."""
    merged_loci = merge_overlapped_loci(loci_data)
    assert merged_loci.shape == (2, 6)
    assert merged_loci["CHR"].tolist() == [1, 2]
    assert merged_loci["START"].tolist() == [100, 100]
    assert merged_loci["END"].tolist() == [600, 600]
    assert merged_loci["LEAD_SNP"].tolist() == ["rs1", "rs6"]
    assert merged_loci["LEAD_SNP_P"].tolist() == [1e-5, 1e-5]
    assert merged_loci["LEAD_SNP_BP"].tolist() == [100, 100]

def test_indep_snps_by_distance(sig_df):
    """Test the indep_snps_by_distance function."""
    indep_snps = indep_snps_by_distance(sig_df)
    assert indep_snps["CHR"].tolist() == [21,22]
    assert indep_snps["BP"].tolist() == [36119111, 18600583]

def test_leadsnp2loci(sig_df):
    """Test the leadsnp2loci function."""
    loci1 = leadsnp2loci(sig_df, if_merge=False)
    loci2 = leadsnp2loci(sig_df, if_merge=True)
    assert ColName.START in loci1.columns
    assert len(loci1) == 9
    assert len(loci2) == 2