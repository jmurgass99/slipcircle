"""Tests para Newmark sliding block."""
import numpy as np
import pytest

from slipcircle.newmark import NewmarkInput, evaluate_newmark


def test_acelerograma_bajo_no_desliza():
    ag = np.full(1000, 0.05)  # 0.05g constante < 0.2g
    r = evaluate_newmark(NewmarkInput(ag, dt_s=0.01, a_y_g=0.2))
    assert r.permanent_displacement_m == 0.0
    assert r.severity == "low"


def test_acelerograma_severo_desliza():
    # 5 segundos de 1g constante con a_y=0.2g
    ag = np.full(500, 1.0)
    r = evaluate_newmark(NewmarkInput(ag, dt_s=0.01, a_y_g=0.2))
    assert r.permanent_displacement_m > 0
    assert r.n_pulses_exceeding > 0


def test_arias_no_negativo():
    rng = np.random.default_rng(42)
    ag = rng.normal(scale=0.3, size=2000)
    r = evaluate_newmark(NewmarkInput(ag, dt_s=0.005, a_y_g=0.15))
    assert r.arias_intensity_m_s > 0


def test_input_invalido():
    with pytest.raises(ValueError):
        evaluate_newmark(NewmarkInput([], dt_s=0.01, a_y_g=0.1))
    with pytest.raises(ValueError):
        evaluate_newmark(NewmarkInput([0.1], dt_s=0, a_y_g=0.1))
