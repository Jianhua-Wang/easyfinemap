"""Tests for the ldref module."""

import os
import pytest


from easyfinemap.ldref import LDRef

PWD = os.path.dirname(os.path.abspath(__file__))

class TestLDRef:
    """Tests for the LDRef class."""

    def test_init(self):
        """Test the initialization of the LDRef class."""
        ldref = LDRef(f"{PWD}/exampledata/EUR.chr21-22", file_type="plink")
        assert ldref.ldref_path == f"{PWD}/exampledata/EUR.chr21-22"
        assert ldref.file_type == "plink"
        assert ldref.temp_dir_path.exists()

    def test_init_error(self):
        """Test the initialization of the LDRef class."""
        with pytest.raises(ValueError):
            LDRef(f"{PWD}/exampledata/EUR.chr21-22", file_type="notype")

    def test_extract(self):
        """Test the extract method of the LDRef class."""
        ldref = LDRef(f"{PWD}/exampledata/EUR.chr21-22", file_type="plink")
        ldref.extract(chrom=21, start=1000000, end=2000000, outprefix="test")
        assert (ldref.temp_dir_path / "test.bim").exists()
        assert (ldref.temp_dir_path / "test.bed").exists()
        assert (ldref.temp_dir_path / "test.fam").exists()