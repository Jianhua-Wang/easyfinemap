"""Tests for the loci module."""

import os
import shutil
from pathlib import Path
import pytest

import pandas as pd

from easyfinemap.loci import Loci
from easyfinemap.constant import ColName

PWD = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()


class TestLoci:
    """Tests for the Loci class."""

    def test_init(self):
        """Test the Loci class."""
        if os.path.exists(f"{CWD}/tmp/loci"):
            shutil.rmtree(f"{CWD}/tmp/loci")
        loci = Loci()
        assert str(loci.tmp_root) == f"{CWD}/tmp/loci"
        for file in Path(f"{PWD}/exampledata/").glob("*loci.txt"):
            if file.is_file():
                os.remove(file)
        for file in Path(f"{PWD}/exampledata/").glob("*leadsnp.txt"):
            if file.is_file():
                os.remove(file)

    def test_identify_indep_loci(self, sumstats_data, clean_ld_panel):
        """Test the identify_indep_loci function."""
        loci = Loci()
        leadsnp, indep_loci = loci.identify_indep_loci(
            sumstats_data,
            method='conditional',
            ldref=f'{clean_ld_panel}.chr{{chrom}}',
            sample_size=1000,
            use_ref_EAF=True,
        )  # type: ignore
        print(leadsnp)
        print(indep_loci)
        assert indep_loci.shape[0] == 4
        assert leadsnp.shape[0] == 4
        loci.identify_indep_loci(
            sumstats_data,
            method='conditional',
            ldref=f'{clean_ld_panel}.chr{{chrom}}',
            sample_size=1000,
            use_ref_EAF=True,
            outprefix=f"{PWD}/exampledata/conditional",
        )
        assert os.path.exists(f"{PWD}/exampledata/conditional.loci.txt")
        assert os.path.exists(f"{PWD}/exampledata/conditional.leadsnp.txt")
        with pytest.raises(ValueError):
            loci.identify_indep_loci(
                sumstats_data,
                method='conditional',
                sample_size=1000,
                use_ref_EAF=True,
                outprefix=f"{PWD}/exampledata/conditional",
            )
        with pytest.raises(ValueError):
            loci.identify_indep_loci(
                sumstats_data,
                method='conditional',
                ldref=f'{clean_ld_panel}.chr{{chrom}}',
                use_ref_EAF=True,
                outprefix=f"{PWD}/exampledata/conditional",
            )
        with pytest.raises(ValueError):
            loci.identify_indep_loci(
                sumstats_data,
                method='nomethod',
                ldref=f'{clean_ld_panel}.chr{{chrom}}',
                sample_size=1000,
                outprefix=f"{PWD}/exampledata/conditional",
            )
        with pytest.raises(ValueError):
            loci.identify_indep_loci(
                sumstats_data,
                method='conditional',
                ldref=f'{clean_ld_panel}.chr{{chrom}}',
                sample_size=1000,
                use_ref_EAF=False,
                outprefix=f"{PWD}/exampledata/conditional",
            )
        loci.identify_indep_loci(
            sumstats_data,
            method='conditional',
            ldref=f'{clean_ld_panel}.chr{{chrom}}',
            sample_size=1000,
            use_ref_EAF=True,
            outprefix=f"{PWD}/exampledata/conditional",
            if_merge=True,
        )

    def test_merge_overlapped_loci(self, loci_data):
        """Test the merge_overlapped_loci function."""
        loci = Loci()
        merged_loci = loci.merge_overlapped_loci(loci_data)
        assert merged_loci.shape == (2, 6)
        assert merged_loci["CHR"].tolist() == [1, 2]
        assert merged_loci["START"].tolist() == [100, 100]
        assert merged_loci["END"].tolist() == [600, 600]
        assert merged_loci["LEAD_SNP"].tolist() == ["rs1", "rs6"]
        assert merged_loci["LEAD_SNP_P"].tolist() == [1e-5, 1e-5]
        assert merged_loci["LEAD_SNP_BP"].tolist() == [100, 100]

    def test_indep_snps_by_distance(self, sumstats_data):
        """Test the indep_snps_by_distance function."""
        loci = Loci()
        loci.identify_indep_loci(sumstats_data, method='distance', outprefix=f"{PWD}/exampledata/distance")
        assert os.path.exists(f"{PWD}/exampledata/distance.loci.txt")
        assert os.path.exists(f"{PWD}/exampledata/distance.leadsnp.txt")

    def test_indep_snps_by_ldclumping(self, sumstats_data, clean_ld_panel):
        """Test the indep_snps_by_ldclumping function."""
        loci = Loci()
        loci.identify_indep_loci(
            sumstats_data,
            ldref=f'{clean_ld_panel}.chr{{chrom}}',
            method='clumping',
            outprefix=f"{PWD}/exampledata/ldclumping",
        )
        assert os.path.exists(f"{PWD}/exampledata/ldclumping.loci.txt")
        assert os.path.exists(f"{PWD}/exampledata/ldclumping.leadsnp.txt")
        with pytest.raises(ValueError):
            loci.identify_indep_loci(
                sumstats_data,
                method='clumping',
                outprefix=f"{PWD}/exampledata/ldclumping",
            )

    def test_leadsnp2loci(self, sig_df):
        """Test the leadsnp2loci function."""
        loci = Loci()
        loci1 = loci.leadsnp2loci(sig_df, if_merge=False)
        loci2 = loci.leadsnp2loci(sig_df, if_merge=True)
        assert ColName.START in loci1.columns
        assert len(loci1) == 9
        assert len(loci2) == 2
