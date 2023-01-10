"""Tools for LD matrix and fine-mapping."""

import logging
import shutil
from typing import Optional


class Tools:
    """check if tools are installed and return their path."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger("Tools")

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
            self.logger.error(f"{tool} is not installed. Please install it first and make sure it is in your PATH.")
            raise ValueError(f"{tool} is not installed. Please install it first and make sure it is in your PATH.")

    @property
    def plink(self):
        """Check if plink is installed."""
        return self._check_tool("plink")

    @property
    def bcftools(self):
        """Check if bcftools is installed."""
        return self._check_tool("bcftools")

    @property
    def gcta(self):
        """Check if gcta is installed."""
        return self._check_tool("gcta64")

    @property
    def finemap(self):
        """Check if finemap is installed."""
        return self._check_tool("finemap")

    @property
    def paintor(self):
        """Check if paintor is installed."""
        return self._check_tool("PAINTOR")

    @property
    def caviarbf(self):
        """Check if caviarbf is installed."""
        return self._check_tool("caviarbf")

    @property
    def model_search(self):
        """Check if model_search is installed."""
        return self._check_tool("model_search")
