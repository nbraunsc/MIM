"""
Unit and regression test for the mim package.
"""

# Import package, test suite, and other packages as needed
import mim
import pytest
import sys

def test_mim_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mim" in sys.modules
