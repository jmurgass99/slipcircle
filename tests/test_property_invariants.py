"""Property-based tests con hypothesis (Mejora C).

Invariantes fisicos que el factor de seguridad debe satisfacer:

1. FS > 0 para todo suelo geometricamente valido.
2. FS sismico (kh > 0) <= FS estatico.
3. FS monotonamente creciente con la cohesion (c' mayor -> FS mayor).
4. FS monotonamente creciente con el angulo de friccion (phi' mayor -> FS mayor).
5. Pendiente mas suave (menor altura para el mismo ancho) -> FS mayor.
6. Bishop y Janbu deben dar resultados del mismo orden de magnitud (factor 2).
"""
from __future__ import annotations
import pytest

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings, strategies as st, HealthCheck

from slopestabilityjm.geometry import Slope, CircularSurface
from slopestabilityjm.soil_profile import SoilLayer
from slopestabilityjm.bishop import BishopAnalysis
from slopestabilityjm.janbu import JanbuAnalysis
from slopestabilityjm.seismic import PseudoStaticAnalysis


def _slope(c: float, phi: float, gamma: float = 18.0,
           height: float = 10.0) -> Slope:
    """Construye 1V:2H talud con altura ``height`` y suelo (c, phi, gamma)."""
    layer = SoilLayer("S", cohesion=c, friction_angle=phi, unit_weight=gamma)
    return Slope(
        profile=[(0.0, 0.0), (2.0 * height, height), (2.0 * height + 10.0, height)],
        layers=[layer], layer_depths=[max(20.0, 3 * height)],
    )


def _surface(height: float = 10.0) -> CircularSurface:
    """Circulo razonable para la geometria de _slope."""
    return CircularSurface(cx=1.5 * height, cy=2.2 * height, R=2.4 * height)


_SETTINGS = settings(
    max_examples=15, deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)


@_SETTINGS
@given(c=st.floats(5.0, 50.0), phi=st.floats(10.0, 35.0))
def test_fs_strictly_positive(c: float, phi: float) -> None:
    """Para suelo no patologico, Bishop debe dar FS finito y positivo."""
    slope = _slope(c, phi)
    r = BishopAnalysis(slope, n_slices=20).factor_of_safety(_surface())
    assert r["FS"] is not None
    assert r["FS"] > 0


@_SETTINGS
@given(c=st.floats(10.0, 40.0), phi=st.floats(15.0, 30.0),
       kh=st.floats(0.05, 0.30))
def test_seismic_fs_le_static_fs(c: float, phi: float, kh: float) -> None:
    """FS sismico (kh > 0) <= FS estatico para mismo talud/superficie."""
    slope = _slope(c, phi)
    surf = _surface()
    fs_static = BishopAnalysis(slope, n_slices=20).factor_of_safety(surf)["FS"]
    fs_seis = PseudoStaticAnalysis(slope, kh=kh).factor_of_safety(surf)["FS"]
    assert fs_seis <= fs_static + 1e-6


@_SETTINGS
@given(c_low=st.floats(5.0, 20.0), c_step=st.floats(5.0, 30.0),
       phi=st.floats(15.0, 30.0))
def test_fs_monotone_in_cohesion(c_low: float, c_step: float, phi: float) -> None:
    """Aumentar c' siempre aumenta FS (todo lo demas constante)."""
    surf = _surface()
    fs_low = BishopAnalysis(_slope(c_low, phi), n_slices=20).factor_of_safety(surf)["FS"]
    fs_hi = BishopAnalysis(_slope(c_low + c_step, phi), n_slices=20).factor_of_safety(surf)["FS"]
    assert fs_hi > fs_low


@_SETTINGS
@given(phi_low=st.floats(10.0, 20.0), phi_step=st.floats(5.0, 15.0),
       c=st.floats(10.0, 30.0))
def test_fs_monotone_in_friction(phi_low: float, phi_step: float, c: float) -> None:
    """Aumentar phi' siempre aumenta FS."""
    surf = _surface()
    fs_low = BishopAnalysis(_slope(c, phi_low), n_slices=20).factor_of_safety(surf)["FS"]
    fs_hi = BishopAnalysis(_slope(c, phi_low + phi_step), n_slices=20).factor_of_safety(surf)["FS"]
    assert fs_hi > fs_low


@_SETTINGS
@given(c=st.floats(15.0, 35.0), phi=st.floats(18.0, 30.0))
def test_bishop_and_janbu_similar_order(c: float, phi: float) -> None:
    """Bishop y Janbu para mismo talud/superficie: ratio en [0.5, 2.0]."""
    slope = _slope(c, phi)
    surf = _surface()
    fs_b = BishopAnalysis(slope, n_slices=25).factor_of_safety(surf)["FS"]
    fs_j = JanbuAnalysis(slope, n_slices=25).factor_of_safety(surf)["FS"]
    ratio = fs_b / fs_j
    assert 0.5 < ratio < 2.0, f"Bishop/Janbu ratio fuera de rango: {ratio}"


@_SETTINGS
@given(kh_low=st.floats(0.02, 0.10), kh_step=st.floats(0.05, 0.20),
       c=st.floats(15.0, 30.0), phi=st.floats(18.0, 28.0))
def test_higher_kh_lowers_fs(kh_low: float, kh_step: float,
                              c: float, phi: float) -> None:
    """kh mayor produce FS sismico menor (monotonia decreciente)."""
    slope = _slope(c, phi)
    surf = _surface()
    fs_low_kh = PseudoStaticAnalysis(slope, kh=kh_low).factor_of_safety(surf)["FS"]
    fs_hi_kh = PseudoStaticAnalysis(slope, kh=kh_low + kh_step).factor_of_safety(surf)["FS"]
    assert fs_hi_kh < fs_low_kh + 1e-6
