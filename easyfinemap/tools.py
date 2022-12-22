"""Tools for LD matrix and fine-mapping."""

import shutil
import sys
from typing import Optional, Union

from easyfinemap.logger import logger


class Tools:
    """check if tools are installed and return their path."""
    def __init__(self, log_level: str = "WARNING"):
        """Initialize."""
        self.logger = logger
        self.logger.setLevel(log_level)

    def _check_tool(self, tool: str) -> Optional[str]:
        """
        Check if tool is installed.

        Parameters
        ----------
        tool : str
            The tool name.

        Returns
        -------
        Optional[str]
            The path of the tool if it is installed, otherwise exit the program.
        """
        tool_path = shutil.which(tool)
        if tool_path:
            return tool_path
        else:
            self.logger.error(f"{tool} is not installed. Please install it first or make sure it is in your PATH.")
            raise ValueError(f"{tool} is not installed. Please install it first or make sure it is in your PATH.")

    @property
    def plink(self):
        """Check if plink is installed."""
        return self._check_tool("plink")

    @property
    def bcftools(self):
        """Check if bcftools is installed."""
        return self._check_tool("bcftools")





