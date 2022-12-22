"""Tests for the loci module."""

import pytest

from easyfinemap.loci import merge_overlapped_loci

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

