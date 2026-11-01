# -*- coding: utf-8 -*-
"""Spencer y Morgenstern-Price: forma cerrada en φ = 0, coherencia con Bishop y entre sí.

Hasta 2026-10-31 los dos métodos (a) usaban en Q_i el numerador de Bishop (c·b y W − u·b)
en lugar del de Spencer (1967) Ec. 7 (c·l y W·cos α − u·l, con l = b·sec α), y (b)
alternaban raíces en FS y en θ/λ con un bracket que incluía los polos de Q_i, de modo
que el resultado dependía del número de dovelas (θ = 0 sin equilibrio de momentos, o
FS = 8,7 con θ = −46°) y la bandera ``converged`` era siempre False. Estos tests anclan:

1. Forma cerrada no drenada (φ' = 0): FS = c·R·L_arco / Σ W·(x_m − c_x). Con φ = 0 la
   ecuación de momentos de Spencer se reduce a Σ(c·l/FS − W·sin α) = 0, que es
   exactamente esa expresión; antes salía Σ c·b (sin sec α) y el FS era menor.
2. Círculo de Duncan y Wright: Spencer y M-P convergen, están a menos del 2 % de Bishop
   y varían menos del 1 % entre 30 y 60 dovelas.
3. M-P con f constante ≡ Spencer (misma θ para todas las dovelas, Q por el punto medio
   de la base).
"""
from __future__ import annotations

from math import asin, radians, tan

import pytest

from slipcircle import (
    BishopAnalysis,
    CircularSurface,
    MorgensternPriceAnalysis,
    Slope,
    SoilLayer,
    SpencerAnalysis,
    constant_f,
)


def _clay_case():
    slope = Slope(profile=[(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)],
                  layers=[SoilLayer("arcilla", cohesion=25.0, friction_angle=0.0, unit_weight=18.0)],
                  layer_depths=[40.0])
    return slope, CircularSurface(15.0, 22.0, 24.0)


def _duncan_case():
    H, beta = 10.0, radians(30.0)
    dx = H / tan(beta)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    return slope, CircularSurface(3.189, 16.909, 18.091)


@pytest.mark.safety_critical
def test_phi_zero_spencer_and_mp_match_undrained_closed_form():
    slope, surf = _clay_case()
    n = 200
    sl = slope.discretize_circle(surf, n)
    x_lo = sl[0]["x_mid"] - sl[0]["width"] / 2
    x_hi = sl[-1]["x_mid"] + sl[-1]["width"] / 2
    L_arc = surf.R * abs(asin((x_hi - surf.cx) / surf.R) - asin((x_lo - surf.cx) / surf.R))
    fs_closed = 25.0 * surf.R * L_arc / sum(s["weight_kN"] * (s["x_mid"] - surf.cx) for s in sl)
    rs = SpencerAnalysis(slope, n_slices=n).factor_of_safety(surf)
    rm = MorgensternPriceAnalysis(slope, n_slices=n).factor_of_safety(surf)
    assert rs["converged"] and rm["converged"]
    assert rs["FS"] == pytest.approx(fs_closed, rel=2e-3)
    assert rm["FS"] == pytest.approx(fs_closed, rel=2e-3)


def test_duncan_wright_spencer_mp_converge_near_bishop_and_stable_in_slices():
    slope, surf = _duncan_case()
    fs_b = BishopAnalysis(slope, n_slices=40).factor_of_safety(surf)["FS"]
    fs_s, fs_m = [], []
    for n in (30, 40, 60):
        rs = SpencerAnalysis(slope, n_slices=n).factor_of_safety(surf)
        rm = MorgensternPriceAnalysis(slope, n_slices=n).factor_of_safety(surf)
        assert rs["converged"], rs
        assert rm["converged"], rm
        assert 10.0 < rs["theta_deg"] < 35.0          # inclinación física de las fuerzas entre dovelas
        fs_s.append(rs["FS"]); fs_m.append(rm["FS"])
    for fs in fs_s + fs_m:
        assert fs == pytest.approx(fs_b, rel=0.02)
    assert max(fs_s) - min(fs_s) < 0.01 * fs_b
    assert max(fs_m) - min(fs_m) < 0.01 * fs_b


def test_mp_constant_f_equals_spencer():
    slope, surf = _duncan_case()
    rs = SpencerAnalysis(slope, n_slices=40).factor_of_safety(surf)
    rm = MorgensternPriceAnalysis(slope, n_slices=40, f=constant_f).factor_of_safety(surf)
    assert rm["FS"] == pytest.approx(rs["FS"], abs=2e-4)
    assert rm["lambda"] == pytest.approx(tan(radians(rs["theta_deg"])), rel=1e-2)


def test_result_exposes_both_equilibrium_factors():
    slope, surf = _duncan_case()
    r = SpencerAnalysis(slope, n_slices=40).factor_of_safety(surf)
    assert r["FS_force_eq"] == pytest.approx(r["FS_moment_eq"], rel=1e-4)


def test_polyline_sampled_from_circle_reproduces_circular_fs():
    """La rama de polilínea debe dar lo mismo que la circular si la polilínea ES el círculo."""
    import numpy as np
    slope, surf = _duncan_case()
    xs = np.linspace(max(surf.cx - surf.R, 0.0), surf.cx + surf.R, 400)
    ys = surf.cy - np.sqrt(np.maximum(surf.R ** 2 - (xs - surf.cx) ** 2, 0.0))
    inside = ys < np.array([slope.ground_y(x) for x in xs])
    poly = list(zip(xs[inside][::20].tolist(), ys[inside][::20].tolist()))
    for cls in (SpencerAnalysis, MorgensternPriceAnalysis):
        rc = cls(slope, n_slices=40).factor_of_safety(surf)
        rp = cls(slope, n_slices=40).factor_of_safety(poly)
        assert rc["converged"] and rp["converged"], (cls.__name__, rc, rp)
        assert rp["FS"] == pytest.approx(rc["FS"], rel=0.01), cls.__name__


def test_seismic_spencer_converges_and_is_bracketed_by_bishop_pseudostatic():
    """Con k_h = 0,10 Spencer converge y queda a menos del 10 % del Bishop pseudo-estatico.

    La diferencia (1,037 frente a 0,956) es de modelo, no de signo: Spencer aplica W_h en el
    centroide de la dovela y PseudoStaticAnalysis en la base (R*k_h*W*cos a) y ademas usa W
    (no W*(1-k_v)) en el momento motriz. Con W_h en la base, Spencer da 1,004 y una
    reimplementacion directa de la formula documentada de Bishop 1,005 (sondeo 2026-10-31).
    """
    from slipcircle import PseudoStaticAnalysis
    slope, surf = _duncan_case()
    fs_ps = PseudoStaticAnalysis(slope=slope, kh=0.10, kv=0.05, n_slices=40).factor_of_safety(surf)["FS"]
    r = SpencerAnalysis(slope, n_slices=40).factor_of_safety(surf, kh=0.10, kv=0.05)
    assert r["converged"], r
    assert r["FS"] == pytest.approx(fs_ps, rel=0.10)
    assert r["FS"] < SpencerAnalysis(slope, n_slices=40).factor_of_safety(surf)["FS"]
