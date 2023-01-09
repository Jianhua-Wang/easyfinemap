"""Tests for EasyFinemap."""

import os
import shutil
from pathlib import Path
import pytest

import pandas as pd

from easyfinemap.loci import Loci
from easyfinemap.constant import ColName
from easyfinemap.easyfinemap import EasyFinemap

PWD = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()


class TestEasyFinemap:
    """Tests for the EasyFinemap class."""

    def test_init(self):
        """Test the EasyFinemap class."""
        if os.path.exists(f"{CWD}/tmp/easyfinemap"):
            shutil.rmtree(f"{CWD}/tmp/easyfinemap")
        easyfinemap = EasyFinemap()
        assert str(easyfinemap.tmp_root) == f"{CWD}/tmp/easyfinemap"
        for file in Path(f"{PWD}/exampledata/").glob("*loci.txt"):
            if file.is_file():
                os.remove(file)
        for file in Path(f"{PWD}/exampledata/").glob("*leadsnp.txt"):
            if file.is_file():
                os.remove(file)
