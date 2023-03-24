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
            f"{PWD}/exampledata/LDREF/EUR.chr21-22",
            f"{PWD}/exampledata/LDREF/EUR.valid",
            "--file-type",
            "plink",
        ],
    )
    assert result.exit_code == 0
    assert "EUR.valid.chr21.bim" in os.listdir(f"{PWD}/exampledata/LDREF")


def test_get_loci():
    """Test the get-loci command."""
    result = runner.invoke(app, ["get-loci", "--help"])
    assert result.exit_code == 0
    assert "get-loci" in result.stdout
    result = runner.invoke(app, ["get-loci", f"{PWD}/exampledata/noEAF_noMAF.txt.gz", f"{PWD}/exampledata/distance"])
    assert result.exit_code == 0
    assert os.path.exists(f"{PWD}/exampledata/distance.loci.txt")
    result = runner.invoke(app, ["get-loci", "None_file", f"{PWD}/exampledata/distance"])
    assert result.exit_code == 1
