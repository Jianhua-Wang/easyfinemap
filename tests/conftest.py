"""Define the test suite for the sumstats module."""

import pytest

import pandas as pd


@pytest.fixture
def loci_data():
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
    sig_df = {
        "CHR": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2],
        "BP": [100, 200, 300, 400, 500, 100, 200, 300, 400, 500],
        "SNP": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9", "rs10"],
        "P": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
    }
    return pd.DataFrame(sig_df)