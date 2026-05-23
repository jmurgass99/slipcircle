"""Tests FS mínimos por código."""
import pytest

from slopestabilityjm.fs_by_code import fs_minimum, is_safe


def test_5_codes():
    for c in ("EU", "SG", "US", "CO", "AU"):
        fs = fs_minimum(c)
        assert fs > 0


def test_us_estatico_1_5():
    assert fs_minimum("US") == 1.5


def test_eu_lrfd_1_0():
    assert fs_minimum("EU") == 1.0


def test_seismic_lower():
    assert fs_minimum("US", "seismic") < fs_minimum("US")


def test_is_safe():
    assert is_safe(1.6, "US")
    assert not is_safe(1.0, "CO")
