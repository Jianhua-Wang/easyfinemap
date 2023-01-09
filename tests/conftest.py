"""Define the test suite for the sumstats module."""

import os

import numpy as np
import pandas as pd
import pytest

PWD = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def sumstats_data():
    """Load the example sumstats data."""
    sumstats_data = pd.read_csv(f"{PWD}/exampledata/noEAF_noMAF.txt.gz", sep="\t")
    return sumstats_data


@pytest.fixture
def loci_data():
    """Load the example loci data."""
    loci_data = {
        "CHR": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
        "START": [100, 200, 300, 400, 500, 100, 200, 300, 400, 500],
        "END": [200, 300, 400, 500, 600, 200, 300, 400, 500, 600],
        "LEAD_SNP": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9", "rs10"],
        "LEAD_SNP_P": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
        "LEAD_SNP_BP": [100, 200, 300, 400, 500, 100, 200, 300, 400, 500],
    }
    return pd.DataFrame(loci_data)


@pytest.fixture
def sig_df():
    """Load the example sig data."""
    sig_df = pd.read_csv(f"{PWD}/exampledata/sig.txt", sep="\t")
    return sig_df


@pytest.fixture
def mock_sumstat():
    """Mock sumstat data for testing."""
    mock_data = pd.DataFrame(
        {
            "CHR": ["chr1", 1, 1, np.nan, 1, 2, 2, 2, 2],
            "BP": [50, 100, 100, 400, None, 100, 200, 300, 400],
            "P": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0, 1e-4, 1e-3, 1e-2],
            "rsID": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9"],
            "EA": ["G", "A", "T", "A", "A", "A", "A", None, "A"],
            "NEA": ["C", "T", "A", "T", "T", "T", "T", "T", "T"],
            "BETA": [0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4],
            "SE": [0.01, 0.02, 0.03, 0.04, 0.05, 0.01, 0.02, 0.03, 0],
            "EAF": [0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4],
        }
    )
    return mock_data


@pytest.fixture
def dirty_ld_panel():
    """Load the example dirty ld panel data."""
    return f"{PWD}/exampledata/LDREF/EUR.chr21-22"


@pytest.fixture
def clean_ld_panel():
    """Load the example clean ld panel data."""
    return f"{PWD}/exampledata/LDREF/EUR.valid"
