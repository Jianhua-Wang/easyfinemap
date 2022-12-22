"""Define the test suite for the sumstats module."""

import pytest
import os

import pandas as pd

PWD = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def sumstats_data():
    sumstats_data = pd.read_csv(f"{PWD}/exampledata/noEAF_noMAF.txt.gz", sep="\t")
    return pd.DataFrame(sumstats_data)


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
    sig_df = pd.read_csv(f"{PWD}/exampledata/sig.txt", sep="\t")
    return sig_df


@pytest.fixture
def mock_sumstat():
    mock_data = pd.DataFrame(
        {
            "CHR": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "BP": [100, 100, 300, 400, 500, 100, 200, 300, 400],
            "P": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e-5, 1e-4, 1e-3, 1e-2],
            "SNP": ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9"],
            "EA": ["A", "T", "G", "A", "A", "A", "A", "A", "A"],
            "NEA": ["T", "A", "C", "T", "T", "T", "T", "T", "T"],
        }
    )
    return mock_data

