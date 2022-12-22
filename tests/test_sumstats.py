"""Tests for the sumstats module."""

import pytest

from easyfinemap.sumstats import SumStats

@pytest.fixture
def sumstats():
    """Return a SumStats object."""
    return SumStats(
        {
            "CHR": [1, 1, 1, 1, 1],
            "BP": [100, 200, 300, 400, 500],
            "rsID": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "EA": ["A", "A", "A", "A", "A"],
            "NEA": ["T", "T", "T", "T", "T"],
            "EAF": [0.1, 0.2, 0.3, 0.4, 0.5],
            "P": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
            "BETA": [0.1, 0.2, 0.3, 0.4, 0.5],
            "SE": [0.1, 0.1, 0.1, 0.1, 0.1],
        }
    )

def test_sumstats(sumstats):
    """Test the SumStats class."""
    assert sumstats["CHR"].dtype == "int64"
    assert sumstats["BP"].dtype == "int64"
    assert sumstats["EAF"].dtype == "float64"
    assert sumstats["P"].dtype == "float64"
    assert sumstats["BETA"].dtype == "float64"
    assert sumstats["SE"].dtype == "float64"