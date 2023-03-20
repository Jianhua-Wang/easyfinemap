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
    SNPID = 'SNPID'  # unique snpid, chr-bp-sorted(EA,NEA)
    COJO_P = "COJO_P"
    COJO_BETA = "COJO_BETA"
    COJO_SE = "COJO_SE"
    LEAD_SNP = "LEAD_SNP"
    LEAD_SNP_P = "LEAD_SNP_P"
    LEAD_SNP_BP = "LEAD_SNP_BP"
    # posterior probability
    PP_FINEMAP = "PP_FINEMAP"
    PP_ABF = "PP_ABF"
    PP_PAINTOR = "PP_PAINTOR"
    PP_CAVIARBF = "PP_CAVIARBF"
    sumstat_cols = ['CHR', 'BP', 'rsID', 'EA', 'NEA', 'P', 'BETA', 'SE', 'EAF', 'MAF']
    loci_cols = ['CHR', 'START', 'END', 'LEAD_SNP', 'LEAD_SNP_P', 'LEAD_SNP_BP']


# only support autosomes
CHROMS = [i for i in range(1, 24)]
