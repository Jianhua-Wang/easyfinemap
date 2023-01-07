"""Tests for the tools module."""


from easyfinemap.tools import Tools


def test_tools():
    """Test the tools class."""
    tools = Tools()
    assert tools.plink
    assert tools.bcftools
    assert tools.gcta
    assert tools.finemap
    assert tools.paintor
    assert tools.caviarbf
    try:
        tools._check_tool("notool")
    except ValueError:
        assert True
