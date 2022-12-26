"""Tests for the cli module."""

from typer.testing import CliRunner

from easyfinemap.cli import app

runner = CliRunner()

def test_main():
    """Test the main entrypoint."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0

def test_validate_ldref():
    """Test the validate_ldref command."""
    result = runner.invoke(app, ["validate-ldref", "--help"])
    assert result.exit_code == 0
    assert "validate-ldref" in result.stdout
    result = runner.invoke(app, ["validate-ldref", "tests/exampledata/LDREF/hapmap3", "tests/exampledata/LDREF/hapmap3.valid"])