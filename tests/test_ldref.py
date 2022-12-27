"""Tests for the ldref module."""

import os
import shutil

import pytest

from easyfinemap.ldref import LDRef

PWD = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()


class TestLDRef:
    """Tests for the LDRef class."""

    def test_init(self):
        """Test the initialization of the LDRef class."""
        if os.path.exists(f"{CWD}/tmp"):
            shutil.rmtree(f"{CWD}/tmp")
        ldref = LDRef()
        assert ldref.temp_dir_path.exists()

    def test_valid(self, dirty_ld_panel, clean_ld_panel):
        """Test the valid method of the LDRef class."""
        if os.path.exists(f"{clean_ld_panel}.chr1.bim"):
            os.remove(f"{clean_ld_panel}.chr1.bim")
        ldref = LDRef()
        ldref.valid(dirty_ld_panel, clean_ld_panel)
        assert os.path.exists(f"{clean_ld_panel}.chr1.bim")
        with pytest.raises(ValueError):
            ldref.valid(dirty_ld_panel, clean_ld_panel, file_type="vcf")
        ldref.valid(f"{clean_ld_panel}.chr{{chrom}}", f"{clean_ld_panel}.copy", file_type="plink")
        assert os.path.exists(f"{clean_ld_panel}.copy.chr1.bim")
        for chrom in range(1, 23):
            os.remove(f"{clean_ld_panel}.copy.chr{chrom}.bim")
            os.remove(f"{clean_ld_panel}.copy.chr{chrom}.bed")
            os.remove(f"{clean_ld_panel}.copy.chr{chrom}.fam")
        with pytest.raises(FileNotFoundError):
            ldref.valid(f"{clean_ld_panel}.123.chr{{chrom}}", f"{clean_ld_panel}.copy", file_type="plink")
        with pytest.raises(FileNotFoundError):
            ldref.valid(f"{dirty_ld_panel}.123", f"{clean_ld_panel}.copy")

    def test_extract(self, clean_ld_panel):
        """Test the extract method of the LDRef class."""
        if os.path.exists(f"{PWD}/exampledata/LDREF/extract.bim"):
            os.remove(f"{PWD}/exampledata/LDREF/extract.bim")
        ldref = LDRef()
        ldref.extract(
            f"{clean_ld_panel}.chr{{chrom}}", f"{PWD}/exampledata/LDREF/extract", chrom=1, start=1000000, end=2000000
        )
        assert os.path.exists(f"{PWD}/exampledata/LDREF/extract.bim")
        with pytest.raises(FileNotFoundError):
            ldref.extract(
                f"{clean_ld_panel}.chr{{chrom}}",
                f"{PWD}/exampledata/LDREF/extract",
                chrom=30,
                start=1000000,
                end=2000000,
            )
        with pytest.raises(RuntimeError):
            ldref.extract(
                f"{clean_ld_panel}.chr{{chrom}}",
                f"{PWD}/exampledata/LDREF/extract",
                chrom=1,
                start=1000000,
                end=2000000,
                mac=10000,
            )

    def test_clean(self):
        """Test the clean method of the LDRef class."""
        if os.path.exists(f"{PWD}/exampledata/LDREF/clean.bim"):
            os.remove(f"{PWD}/exampledata/LDREF/clean.bim")
        ldref = LDRef()
        ldref.clean(f"{PWD}/exampledata/LDREF/extract", f"{PWD}/exampledata/LDREF/clean", mac=10)
        assert os.path.exists(f"{PWD}/exampledata/LDREF/clean.bim")
        with pytest.raises(ValueError):
            ldref.clean(f"{PWD}/exampledata/LDREF/extract", f"{PWD}/exampledata/LDREF/clean", mac=0)
        with pytest.raises(RuntimeError):
            ldref.clean(f"{PWD}/exampledata/LDREF/extract", f"{PWD}/exampledata/LDREF/clean", mac=10000)
        ldref.clean(f"{PWD}/exampledata/LDREF/extract")
        assert os.path.exists(f"{ldref.temp_dir_path}/extract.clean.bim")
