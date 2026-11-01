from __future__ import annotations

import numpy as np
import pytest

from slipcircle.geometry import CircularSurface, Slope
from slipcircle.janbu_generalized import JanbuGeneralizedAnalysis
from slipcircle.morgenstern_price import MorgensternPriceAnalysis
from slipcircle.reinforcement import GeosyntheticLayer, SoilNail
from slipcircle.search import CriticalSurfaceSearch, grid_search, optimize_circle
from slipcircle.soil_profile import SoilLayer
from slipcircle.spencer import SpencerAnalysis


def test_multi_layer_depth_from_local_surface():
    upper = SoilLayer("Fill", cohesion=5.0, friction_angle=28.0, unit_weight=17.0)
    lower = SoilLayer("StiffClay", cohesion=45.0, friction_angle=22.0, unit_weight=20.0)
    slope = Slope(
        profile=[(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)],
        layers=[upper, lower],
        layer_depths=[3.0, 25.0],
    )
    # At x=25 (crest y=10), y=8.5 is 1.5 m below ground -> upper layer
    assert slope.layer_at_depth(8.5, x=25.0).name == "Fill"
    # At x=25 (crest y=10), y=4.0 is 6.0 m below ground -> lower layer
    assert slope.layer_at_depth(4.0, x=25.0).name == "StiffClay"


def test_janbu_generalized_regularized_thrust_and_polyline_reinforcement():
    layer = SoilLayer("Soil", cohesion=15.0, friction_angle=22.0, unit_weight=18.5)
    slope = Slope(profile=[(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)], layers=[layer], layer_depths=[30.0])
    poly = [(2.0, 0.5), (10.0, -1.5), (18.0, 1.0), (25.0, 10.0)]
    geogrid = GeosyntheticLayer(y=4.0, x_left=0.0, x_right=28.0, T_allow=90.0, tau_bond=35.0)

    r0 = JanbuGeneralizedAnalysis(slope, poly, n_slices=25).compute_fs()
    r1 = JanbuGeneralizedAnalysis(slope, poly, n_slices=25, reinforcements=[geogrid]).compute_fs()
    assert r0["FS"] is not None and r1["FS"] is not None
    assert r1["FS"] > r0["FS"]


def test_spencer_and_mp_polyline_reinforcement():
    layer = SoilLayer("Soil", cohesion=18.0, friction_angle=20.0, unit_weight=18.0)
    slope = Slope(profile=[(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)], layers=[layer], layer_depths=[30.0])
    cx, cy, R = 15.0, 20.0, 22.0
    xs = np.linspace(2.0, 26.0, 18)
    ys = cy - np.sqrt(np.maximum(R * R - (xs - cx) ** 2, 0.0))
    poly = list(zip(xs.tolist(), ys.tolist()))

    geogrid = GeosyntheticLayer(y=4.0, x_left=0.0, x_right=30.0, T_allow=100.0, tau_bond=40.0)
    nail = SoilNail(x_head=8.0, y_head=4.0, inclination_deg=15.0, length=18.0, T_allow=120.0)

    for cls in (SpencerAnalysis, MorgensternPriceAnalysis):
        r_unreinf = cls(slope, n_slices=25).factor_of_safety(poly)
        r_reinf = cls(slope, n_slices=25, reinforcements=[geogrid, nail]).factor_of_safety(poly)
        assert r_reinf["FS"] > r_unreinf["FS"]
        assert r_reinf["M_reinforcement_kNm"] > 0.0


def test_critical_surface_search_methods_and_optimize():
    layer = SoilLayer("Soil", cohesion=20.0, friction_angle=20.0, unit_weight=18.0)
    slope = Slope(profile=[(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)], layers=[layer], layer_depths=[30.0])

    r_sp = grid_search(slope, method="spencer", n_cx=3, n_cy=3, n_R=3, n_slices=16)
    assert r_sp["critical_FS"] is not None and r_sp["critical_FS"] > 0.0
    assert "Spencer" in r_sp["method"]

    r_opt = optimize_circle(
        slope,
        method="bishop",
        n_slices=20,
        initial_surface=CircularSurface(15.0, 20.0, 22.0),
    )
    assert r_opt["critical_FS"] is not None and r_opt["critical_FS"] > 0.0
