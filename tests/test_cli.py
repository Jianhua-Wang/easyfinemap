"""Tests for the cli module."""

import os

import pytest
from typer.testing import CliRunner

from easyfinemap.cli import app

runner = CliRunner()
PWD = os.path.dirname(os.path.abspath(__file__))


def test_main():
    """Test the main entrypoint."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0


def test_validate_ldref():
    """Test the validate_ldref command."""
    result = runner.invoke(app, ["validate-ldref", "--help"])
    assert result.exit_code == 0
    assert "validate-ldref" in result.stdout
    result = runner.invoke(
        app,
        [
            "validate-ldref",
            f"{PWD}/exampledata/LDREF/hapmap3",
            f"{PWD}/exampledata/LDREF/hapmap3.valid",
            "--file-type",
            "plink",
        ],
    )
    assert result.exit_code == 0
    assert "hapmap3.valid.chr1.bim" in os.listdir(f"{PWD}/exampledata/LDREF")


def test_validate_sumstats():
    """Test the validate_sumstats command."""
    result = runner.invoke(app, ["validate-sumstats", "--help"])
    assert result.exit_code == 0
    assert "validate-sumstats" in result.stdout
    result = runner.invoke(
        app, ["validate-sumstats", f"{PWD}/exampledata/noEAF_noMAF.txt.gz", f"{PWD}/exampledata/sumstats.valid.txt"],
    )
    assert result.exit_code == 0
    assert "sumstats.valid.txt" in os.listdir(f"{PWD}/exampledata/")
    os.remove(f"{PWD}/exampledata/sumstats.valid.txt")
    result = runner.invoke(
        app, ["validate-sumstats", f"{PWD}/exampledata/noEAF_noMAF.txt.gz", f"{PWD}/exampledata/sumstats.valid.txt.gz"],
    )
    assert result.exit_code == 0
    assert "sumstats.valid.txt.gz" in os.listdir(f"{PWD}/exampledata/")
    os.remove(f"{PWD}/exampledata/sumstats.valid.txt.gz")
    # with pytest.raises(FileNotFoundError):
    #     result = runner.invoke(
    #         app,
    #         [
    #             "validate-sumstats",
    #             f"{PWD}/exampledata/noEAF_noMAF.1.txt.gz",
    #             f"{PWD}/exampledata/sumstats.valid.txt.gz",
    #         ],
    #     )
