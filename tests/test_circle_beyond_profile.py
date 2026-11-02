"""Círculos que salen más allá de los extremos del perfil: el terreno sigue en horizontal.

Antes la discretización recortaba el círculo al tramo [x_pie, x_fin] del perfil y perdía la
cuña pasiva de un círculo profundo que sale por delante del pie (FS = 0,42 frente a 1,11 a mano).
Con φ = 0, Bishop es exacto en momentos: FS = c_u·L·R / Σ W·(x − x_c).
"""
from __future__ import annotations

import numpy as np
import pytest

from slipcircle import BishopAnalysis, CircularSurface, CriticalSurfaceSearch, Slope, SoilLayer

PROFILE = [(0.0, 0.0), (9.0, 6.0), (29.0, 6.0)]


def _slope() -> Slope:
    return Slope(profile=PROFILE, layers=[SoilLayer("arcilla", cohesion=20.0, friction_angle=0.0,
                                                     unit_weight=20.0)], layer_depths=[30.0])


def _fs_phi0_by_hand(cx: float, cy: float, R: float, cu: float = 20.0, gamma: float = 20.0) -> float:
    xs = np.linspace(cx - R + 1e-9, cx + R - 1e-9, 200_001)
    yb = cy - np.sqrt(R ** 2 - (xs - cx) ** 2)
    yt = np.interp(xs, [p[0] for p in PROFILE], [p[1] for p in PROFILE])
    h = np.clip(yt - yb, 0.0, None)
    m = h > 0
    dx = xs[1] - xs[0]
    overturning = np.sum(gamma * h[m] * dx * (xs[m] - cx))
    arc = np.sum(np.sqrt(1.0 + np.gradient(yb, dx)[m] ** 2) * dx)
    return cu * arc * R / overturning


@pytest.mark.parametrize("circle", [(-0.308, 9.164, 15.0), (5.077, 9.164, 10.855), (6.0, 12.0, 13.0)])
def test_bishop_phi0_matches_hand_moment_balance_including_ground_beyond_toe(circle):
    fs = BishopAnalysis(_slope(), n_slices=60).factor_of_safety(CircularSurface(*circle))["FS"]
    assert fs == pytest.approx(_fs_phi0_by_hand(*circle), rel=0.02)


def test_circle_ground_extent_reaches_beyond_profile_start():
    x_lo, x_hi = _slope().circle_ground_extent(CircularSurface(-0.308, 9.164, 15.0))
    assert x_lo < -12.0 and x_hi > 14.0


def test_deep_search_no_longer_halves_fs_against_taylor():
    """φ = 0, β = 33,7°, N = c_u/(γH) = 0,167; Taylor (D grande) N_s ≈ 0,18 → FS ≈ 0,92."""
    res = CriticalSurfaceSearch(_slope(), n_cx=14, n_cy=12, n_R=12, n_slices=30,
                                require_toe_exit=False, max_depth_below_toe=1.0).search()
    assert 0.85 < res["critical_FS"] < 1.05
