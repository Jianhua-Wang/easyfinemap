"""Tests for the tools module."""


from easyfinemap.tools import Tools


def test_tools():
    """Test the tools class."""
    tools = Tools()
    assert tools.plink
    assert tools.bcftools
    try:
        tools._check_tool("notool")
    except ValueError:
        assert True
