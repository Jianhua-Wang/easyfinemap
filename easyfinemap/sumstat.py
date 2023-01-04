"""Standarize summary statistics for use with EasyFINEMAP."""

import logging

import pandas as pd

from easyfinemap.constant import CHROMS, ColName
from easyfinemap.utils import make_SNPID_unique

pd.options.mode.chained_assignment = None  # type: ignore


class SumStat(pd.DataFrame):
    """extension of pd.Dataframe for standarize summary statistics."""

    def __init__(self, *args, **kwargs):
        """Initialize the SumstatAccessor class."""
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger("Sumstat")

    @property
    def _constructor(self):
        return SumStat

    @property
    def _constructor_sliced(self):
        return pd.Series

    def drop_allna_cols(self):
        """Drop all columns with all missing values."""
        for col in self.columns:
            if self[col].isnull().all():
                del self[col]
        return self

    def _validate(self):
        """Check if the required columns are in the summary statistics."""
        required_cols = [
            ColName.CHR,
            ColName.BP,
            ColName.EA,
            ColName.NEA,
            ColName.P,
            ColName.BETA,
            ColName.SE,
        ]
        missing_cols = [col for col in required_cols if col not in self.columns]
        if len(missing_cols) > 0:
            raise ValueError(f"Missing required columns: {missing_cols}")

    def check_chr(self):
        """Check if chromosome is a valid integer."""
        # remove chr from chromosome column, if exists
        self = self[self[ColName.CHR].notnull()].copy()
        self[ColName.CHR] = self[ColName.CHR].astype("str").str.replace("chr", "")
        # replace X, with 23
        self[ColName.CHR] = self[ColName.CHR].replace("X", 23)
        # turn chromosome column into integer
        self[ColName.CHR] = pd.to_numeric(self[ColName.CHR], errors="coerce")
        self = self[self[ColName.CHR].notnull()]
        self[ColName.CHR] = self[ColName.CHR].astype(int)
        self = self[self[ColName.CHR].isin(CHROMS)]
        return self

    def check_bp(self):
        """Check if bp is a valid integer."""
        # raise NotImplementedError
        self[ColName.BP] = pd.to_numeric(self[ColName.BP], errors="coerce")
        self = self[self[ColName.BP].notnull()]
        self[ColName.BP] = self[ColName.BP].astype(int)
        return self

    def check_ea_nea(self):
        """Check if ea and nea are valid."""
        self[ColName.EA] = self[ColName.EA].str.upper()
        self[ColName.NEA] = self[ColName.NEA].str.upper()
        # remove rows with missing ea or nea
        self = self[self[ColName.EA].notnull() & self[ColName.NEA].notnull()]
        return self

    def check_p(self):
        """Check if p is valid."""
        self[ColName.P] = pd.to_numeric(self[ColName.P], errors="coerce")
        self = self[self[ColName.P].notnull()]
        self[ColName.P] = self[ColName.P].astype(float)
        self = self[(self[ColName.P] > 0) & (self[ColName.P] < 1)]
        return self

    def check_beta(self):
        """Check if beta is valid."""
        self[ColName.BETA] = pd.to_numeric(self[ColName.BETA], errors="coerce")
        self = self[self[ColName.BETA].notnull()]
        self[ColName.BETA] = self[ColName.BETA].astype(float)
        return self

    def check_se(self):
        """Check if se is valid."""
        self[ColName.SE] = pd.to_numeric(self[ColName.SE], errors="coerce")
        self = self[self[ColName.SE].notnull()]
        self[ColName.SE] = self[ColName.SE].astype(float)
        self = self[self[ColName.SE] > 0]
        return self

    def check_snpid(self):
        """Check if snpid is valid."""
        if ColName.SNPID in self.columns:
            self.logger.warning("rewriting SNPID column")
            del self[ColName.SNPID]
        result = self.copy()
        allele_df = make_SNPID_unique(result, remove_duplicates=False)
        result.insert(loc=0, column=ColName.SNPID, value=allele_df[ColName.SNPID].values)  # type: ignore
        result.sort_values(ColName.P, inplace=True)
        result.drop_duplicates(subset=[ColName.SNPID], keep="first", inplace=True)
        result.sort_values([ColName.CHR, ColName.BP], inplace=True)
        result.reset_index(drop=True, inplace=True)
        return result

    def check_eaf(self):
        """Check if eaf is valid."""
        if ColName.EAF in self.columns:
            self[ColName.EAF] = pd.to_numeric(self[ColName.EAF], errors="coerce")
            self = self[self[ColName.EAF].notnull()]
            self[ColName.EAF] = self[ColName.EAF].astype(float)
            self = self[(self[ColName.EAF] > 0) & (self[ColName.EAF] < 1)]
        else:
            pass
        return self

    def check_maf(self):
        """Check if maf is valid."""
        if ColName.EAF in self.columns:
            self[ColName.MAF] = self[ColName.EAF]
        else:
            pass
        if ColName.MAF in self.columns:
            self[ColName.MAF] = pd.to_numeric(self[ColName.MAF], errors="coerce")
            self = self[self[ColName.MAF].notnull()]
            self[ColName.MAF] = self[ColName.MAF].astype(float)
            self[ColName.MAF] = self[ColName.MAF].apply(lambda x: min(x, 1 - x))
        else:
            pass
        return self

    def standarize(self):
        """Standarize the data."""
        self._validate()
        self = self.drop_allna_cols()
        self = self.check_chr()
        self = self.check_bp()
        self = self.check_ea_nea()
        self = self.check_p()
        self = self.check_beta()
        self = self.check_se()
        self = self.check_snpid()
        self = self.check_eaf()
        self = self.check_maf()
        self[ColName.Z] = self[ColName.BETA] / self[ColName.SE]
        self.reset_index(drop=True, inplace=True)
        return self
