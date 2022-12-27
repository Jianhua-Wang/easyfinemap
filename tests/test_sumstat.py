"""Tests for the sumstat module."""

import os
import pytest

from easyfinemap.sumstat import SumStat
from easyfinemap.constant import ColName


class TestSumStat:
    """Test the SumStat class."""

    def test_init(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        assert len(sumstat) == len(mock_sumstat)

    def test_validate(self, mock_sumstat):
        """Test the SumStat class."""
        mock_sumstat = mock_sumstat.drop(ColName.CHR, axis=1)
        sumstat = SumStat(mock_sumstat)
        with pytest.raises(ValueError):
            sumstat._validate()

    def test_drop_allna_cols(self, mock_sumstat):
        """Test the SumStat class."""
        mock_sumstat[ColName.CHR] = None
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.drop_allna_cols()
        assert ColName.CHR not in sumstat.columns

    def test_check_chr(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_chr()
        assert sumstat[ColName.CHR].dtype == "int64"
        assert sumstat[ColName.CHR].min() == 1

    def test_check_bp(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_bp()
        assert sumstat[ColName.BP].dtype == "int64"
        assert len(sumstat) == 8

    def test_check_ea_nea(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_ea_nea()
        assert sumstat[ColName.EA].dtype == "object"
        assert sumstat[ColName.NEA].dtype == "object"
        assert len(sumstat) == 8

    def test_check_snpid(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_ea_nea().check_snpid()
        assert sumstat[ColName.SNPID].dtype == "object"
        mock_sumstat[ColName.SNPID] = mock_sumstat[ColName.RSID]
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_ea_nea().check_snpid()
        assert sumstat[ColName.SNPID].dtype == "object"

    def test_check_p(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_p()
        assert sumstat[ColName.P].dtype == "float64"
        assert len(sumstat) == 8

    def test_check_beta(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_beta()
        assert sumstat[ColName.BETA].dtype == "float64"

    def test_check_se(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_se()
        assert sumstat[ColName.SE].dtype == "float64"

    def test_check_eaf(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_eaf()
        assert sumstat[ColName.EAF].dtype == "float64"
        del mock_sumstat[ColName.EAF]
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_eaf()
        assert ColName.EAF not in sumstat.columns

    def test_check_maf(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_maf()
        assert sumstat[ColName.MAF].dtype == "float64"
        mock_sumstat[ColName.MAF] = mock_sumstat[ColName.EAF]
        del mock_sumstat[ColName.EAF]
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_maf()
        assert sumstat[ColName.MAF].dtype == "float64"
        del mock_sumstat[ColName.MAF]
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.check_maf()
        assert ColName.MAF not in sumstat.columns

    def test_standarize(self, mock_sumstat):
        """Test the SumStat class."""
        sumstat = SumStat(mock_sumstat)
        sumstat = sumstat.standarize()
        assert sumstat[ColName.BETA].dtype == "float64"
        assert sumstat[ColName.SE].dtype == "float64"
        assert len(sumstat) == 3
