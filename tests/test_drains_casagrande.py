"""Tests para drenes Casagrande."""
import pytest
from slopestabilityjm.drains_casagrande import (
    DrainDesignInput, design_drains, time_to_drawdown_s,
)


def _typical_input(**overrides):
    base = dict(k_m_s=1e-5, slope_height_m=15.0, target_drawdown_m=5.0,
                drain_length_m=20.0, drain_spacing_m=5.0)
    base.update(overrides)
    return DrainDesignInput(**base)


def test_basic_design():
    r = design_drains(_typical_input())
    assert r.n_drains_required >= 1
    assert r.Q_per_drain_m3_s > 0
    assert r.total_Q_m3_day > 0


def test_higher_k_more_flow():
    a = design_drains(_typical_input(k_m_s=1e-5))
    b = design_drains(_typical_input(k_m_s=1e-4))
    assert b.Q_per_drain_m3_s > a.Q_per_drain_m3_s


def test_multi_code_accepted():
    for code in ("EU", "USA", "NSR", "AS", "SG"):
        r = design_drains(_typical_input(code=code))
        assert r.n_drains_required >= 1
    with pytest.raises(ValueError):
        design_drains(_typical_input(code="ZZ"))


def test_time_to_drawdown():
    t = time_to_drawdown_s(k_m_s=1e-5, drawdown_m=2.0, distance_m=10.0)
    assert t > 0


def test_invalid_inputs():
    with pytest.raises(ValueError):
        design_drains(_typical_input(k_m_s=0))
    with pytest.raises(ValueError):
        design_drains(_typical_input(target_drawdown_m=0))
    with pytest.raises(ValueError):
        time_to_drawdown_s(0, 1, 1)
