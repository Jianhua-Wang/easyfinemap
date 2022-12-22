"""Defines the SumStats class, which is used to store summary statistics."""

import sys
import pandas as pd

from easyfinemap.logger import logger
from easyfinemap.constant import _required_cols


class SumStats(pd.DataFrame):
    """Summary statistics class."""

    _required_cols = _required_cols

    def __init__(self, *args, **kwargs):
        """Initialize."""
        super().__init__(*args, **kwargs)
        self.logger = logger
        self.logger.setLevel("INFO")

    def _check_columns(self):
        """Check if the required columns are in the summary statistics."""
        for col in self._required_cols.keys():
            if col not in self.columns:
                self.logger.error(f"{col} is not in the summary statistics.")
                sys.exit(1)

    def _check_type(self):
        """Check if the column types are correct."""
        for col in self._required_cols.keys():
            dtype = self._required_cols[col]["dtype"]
            if self[col].dtype != dtype:
                self.logger.warning(f"{col} should be in {dtype} format.")
                self.logger.warning(f"Converting {col} to {dtype} format.")
                if self._column_types[col] in self._numeric_columns:
                    self[col] = pd.to_numeric(self[col], errors='coerce')
                else:
                    self[col] = self[col].astype(dtype, errors='ignore')

    def _remove_nan(self):
        """Remove rows with NaN values."""
        for col in self._not_allow_nan_columns:
            na_count = self[col].isnull().sum()
            if na_count:
                self.logger.warning(f"{col} contains {na_count} NaN values.")
                self.logger.warning(f"Removing rows with NaN values in {col}.")
                self.dropna(subset=[col], inplace=True)

    def _replace_chr_x(self):
        """Replace "X" with 23 in the CHR column."""
        self.loc[self["CHR"].isin(["x", "X"]), "CHR"] = 23

    def _check_chr(self):
        """Check if the CHR column is in the correct format."""
        # if chr column starts with "chr", remove "chr"
        self.logger.info("Checking CHR column.")
        if self["CHR"].dtype != "int64" and self["CHR"].str.startswith("chr").any():
            self["CHR"] = self["CHR"].astype('string').str.replace("chr", "", regex=False)
        self._replace_chr_x()

    def validate(self):
        """
        Validate the summary statistics.
        1. Check if the required columns are in the summary statistics.
        2. Check if the column types are correct.
        """
        self._check_columns()
        self._check_type()
        self._check_chr()
        self._remove_nan()

    @property
    def _constructor(self):
        return SumStats

    # @property
    # def _constructor_sliced(self):
    #     return SumStats