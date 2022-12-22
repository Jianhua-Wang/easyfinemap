"""Define constants used in the package."""

class ColName:
    """Define column names."""

    CHR = "CHR"
    BP = "BP"
    RSID = "rsID"
    EA = "EA"
    NEA = "NEA"
    P = "P"
    BETA = "BETA"
    SE = "SE"
    EAF = "EAF"
    MAF = "MAF"
    N = "N"
    Z = "Z"
    INFO = "INFO"
    START = "START"
    END = "END"
    SNPID = 'SNPID'
    LEAD_SNP = "LEAD_SNP"
    LEAD_SNP_P = "LEAD_SNP_P"
    LEAD_SNP_BETA = "LEAD_SNP_BETA"
    LEAD_SNP_SE = "LEAD_SNP_SE"
    LEAD_SNP_EAF = "LEAD_SNP_EAF"
    LEAD_SNP_MAF = "LEAD_SNP_MAF"
    LEAD_SNP_N = "LEAD_SNP_N"
    LEAD_SNP_Z = "LEAD_SNP_Z"
    LEAD_SNP_INFO = "LEAD_SNP_INFO"
    LEAD_SNP_CHR = "LEAD_SNP_CHR"
    LEAD_SNP_BP = "LEAD_SNP_BP"
    LEAD_SNP_RSID = "LEAD_SNP_RSID"
    LEAD_SNP_EA = "LEAD_SNP_EA"
    LEAD_SNP_NEA = "LEAD_SNP_NEA"
    sumstat_cols = ['CHR', 'BP', 'rsID', 'EA', 'NEA', 'P', 'BETA', 'SE', 'EAF', 'MAF']
    loci_cols = ['CHR', 'START', 'END', 'LEAD_SNP', 'LEAD_SNP_P', 'LEAD_SNP_BP']

_required_cols = {
    "CHR": {
        "name": "CHR",
        "description": "Chromosome",
        "dtype": "int64",
        "max": 23,
        "min": 1,
        "allow_nan": False,
        "nonzero": True,
    },
    "BP": {
        "name": "BP",
        "description": "Base pair position",
        "dtype": "int64",
        "max": None,
        "min": 1,
        "allow_nan": False,
        "nonzero": True,
    },
    "rsID": {
        "name": "rsID",
        "description": "SNP ID",
        "dtype": "string",
        "max": None,
        "min": None,
        "allow_nan": True,
        "nonzero": False,
    },
    "EA": {
        "name": "EA",
        "description": "Effect allele",
        "dtype": "string",
        "max": None,
        "min": None,
        "allow_nan": False,
        "nonzero": False,
    },
    "NEA": {
        "name": "NEA",
        "description": "Non-effect allele",
        "dtype": "string",
        "max": None,
        "min": None,
        "allow_nan": False,
        "nonzero": False,
    },
    "P": {
        "name": "P",
        "description": "P-value",
        "dtype": "float64",
        "max": 1,
        "min": 0,
        "allow_nan": False,
        "nonzero": False,
    },
    "BETA": {
        "name": "BETA",
        "description": "Effect size",
        "dtype": "float64",
        "max": None,
        "min": None,
        "allow_nan": False,
        "nonzero": False,
    },
    "SE": {
        "name": "SE",
        "description": "Standard error",
        "dtype": "float64",
        "max": None,
        "min": 0,
        "allow_nan": False,
        "nonzero": True,
    },
}

_optional_cols = {
    "EAF": {
        "name": "EAF",
        "description": "Effect allele frequency",
        "dtype": "float64",
        "max": 1,
        "min": 0,
        "allow_nan": True,
        "nonzero": False,
    },
    "MAF": {
        "name": "MAF",
        "description": "Minor allele frequency",
        "dtype": "float64",
        "max": 1,
        "min": 0,
        "allow_nan": True,
        "nonzero": False,
    },
    "N": {
        "name": "N",
        "description": "Sample size",
        "dtype": "int64",
        "max": None,
        "min": 1,
        "allow_nan": True,
        "nonzero": True,
    },
    "Z": {
        "name": "Z",
        "description": "Z-score",
        "dtype": "float64",
        "max": None,
        "min": None,
        "allow_nan": True,
        "nonzero": False,
    },
    "INFO": {
        "name": "INFO",
        "description": "Imputation quality score",
        "dtype": "float64",
        "max": 1,
        "min": 0,
        "allow_nan": True,
        "nonzero": False,
    },
}

