"""Tests Bray-Travasarou."""
import pytest

from slopestabilityjm.bray_travasarou import (
    BrayTravasarouInput, evaluate_bray_travasarou,
)


def test_input_invalido():
    with pytest.raises(ValueError):
        evaluate_bray_travasarou(BrayTravasarouInput(
            k_y=0, Sa_1_5Ts_g=0.3, Mw=7.0,
        ))


def test_ky_alto_reduce_desplazamiento():
    low = evaluate_bray_travasarou(BrayTravasarouInput(
        k_y=0.05, Sa_1_5Ts_g=0.4, Mw=7.0, Ts_seconds=0.5,
    ))
    high = evaluate_bray_travasarou(BrayTravasarouInput(
        k_y=0.30, Sa_1_5Ts_g=0.4, Mw=7.0, Ts_seconds=0.5,
    ))
    assert high.median_displacement_cm < low.median_displacement_cm


def test_probabilidad_consistente():
    r = evaluate_bray_travasarou(BrayTravasarouInput(
        k_y=0.10, Sa_1_5Ts_g=0.5, Mw=7.5, Ts_seconds=0.5,
    ))
    assert 0 <= r.P_D_gt_30cm <= r.P_D_gt_5cm <= 1.0
