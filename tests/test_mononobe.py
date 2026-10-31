"""Tests Mononobe pseudo-estático."""
from slipcircle.mononobe import MononobeInput, evaluate_mononobe


def test_FS_seismic_menor_static():
    r = evaluate_mononobe(MononobeInput(
        slope_height_m=10, slope_angle_deg=30,
        gamma_kN_m3=18, c_kPa=15, phi_deg=30,
        kh=0.15,
    ))
    assert r.FS_seismic < r.FS_static


def test_sin_sismo_iguales():
    r = evaluate_mononobe(MononobeInput(
        slope_height_m=10, slope_angle_deg=20,
        gamma_kN_m3=18, c_kPa=10, phi_deg=30,
        kh=0.0,
    ))
    import pytest
    assert r.FS_seismic == pytest.approx(r.FS_static, rel=1e-3)


def test_psi_calculado():
    r = evaluate_mononobe(MononobeInput(
        slope_height_m=10, slope_angle_deg=20,
        gamma_kN_m3=18, c_kPa=10, phi_deg=30,
        kh=0.2, kv=0.05,
    ))
    assert r.psi_deg > 0
