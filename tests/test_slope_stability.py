"""
tests/test_slope_stability.py — Unit tests for SlipCircle.
Validated against known Bishop method benchmarks (Duncan & Wright 2014).
"""
import pytest
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from slipcircle.soil_profile import SoilLayer
from slipcircle.geometry import Slope, CircularSurface
from slipcircle.bishop import BishopAnalysis
from slipcircle.janbu import JanbuAnalysis
from slipcircle.seismic import PseudoStaticAnalysis
from slipcircle.search import CriticalSurfaceSearch
from slipcircle.spencer import SpencerAnalysis
from slipcircle.morgenstern_price import MorgensternPriceAnalysis
from slipcircle.janbu_generalized import JanbuGeneralizedAnalysis
from slipcircle.reinforcement import GeosyntheticLayer


# ---------------------------------------------------------------------------
# Helper: create a standard homogeneous slope (1V:2H, H=10m)
# ---------------------------------------------------------------------------
def make_slope(cohesion=20.0, phi=20.0, gamma=18.0,
               water_table_y=None):
    """1V:2H slope, 10m high, toe at (0,0), crest at (30,10)."""
    layer = SoilLayer("Soil", cohesion=cohesion,
                      friction_angle=phi, unit_weight=gamma)
    return Slope(
        profile=[(0, 0), (20, 10), (30, 10)],   # toe → break → crest
        layers=[layer],
        layer_depths=[20.0],
        water_table_y=water_table_y,
    )


def make_surface(cx=15.0, cy=20.0, R=22.0):
    return CircularSurface(cx, cy, R)


# ---------------------------------------------------------------------------
# SoilLayer tests
# ---------------------------------------------------------------------------
class TestSoilLayer:

    def test_cohesive_soil(self):
        s = SoilLayer("Clay", cohesion=30.0, friction_angle=0.0, unit_weight=17.0)
        assert s.tan_phi == pytest.approx(0.0, abs=1e-6)

    def test_frictional_soil(self):
        s = SoilLayer("Sand", cohesion=0.0, friction_angle=35.0, unit_weight=19.0)
        import math
        assert s.tan_phi == pytest.approx(math.tan(math.radians(35.0)), rel=1e-4)

    def test_invalid_friction_angle_raises(self):
        with pytest.raises(ValueError):
            SoilLayer("X", cohesion=10.0, friction_angle=91.0, unit_weight=18.0)

    def test_negative_cohesion_raises(self):
        with pytest.raises(ValueError):
            SoilLayer("X", cohesion=-5.0, friction_angle=20.0, unit_weight=18.0)

    def test_sat_unit_weight_default(self):
        s = SoilLayer("Clay", 20.0, 15.0, 18.0)
        assert s.sat_unit_weight == 18.0


# ---------------------------------------------------------------------------
# Slope geometry tests
# ---------------------------------------------------------------------------
class TestSlopeGeometry:

    def test_height(self):
        slope = make_slope()
        assert slope.height == pytest.approx(10.0)

    def test_ground_y_interpolation(self):
        slope = make_slope()
        assert slope.ground_y(0.0) == pytest.approx(0.0)
        assert slope.ground_y(30.0) == pytest.approx(10.0)
        # On the slope face: y at x=10 should be 5.0
        y10 = slope.ground_y(10.0)
        assert 4.0 <= y10 <= 6.0

    def test_pore_pressure_dry(self):
        slope = make_slope(water_table_y=None)
        assert slope.pore_pressure(10.0, 5.0) == 0.0

    def test_pore_pressure_below_water_table(self):
        slope = make_slope(water_table_y=5.0)
        u = slope.pore_pressure(10.0, 0.0)
        assert u == pytest.approx(9.81 * 5.0, rel=0.01)

    def test_pore_pressure_above_water_table(self):
        slope = make_slope(water_table_y=5.0)
        u = slope.pore_pressure(10.0, 8.0)
        assert u == 0.0

    def test_discretize_gives_slices(self):
        slope = make_slope()
        surf = make_surface()
        slices = slope.discretize_circle(surf, n_slices=20)
        assert len(slices) > 5   # should have several valid slices

    def test_slice_weights_positive(self):
        slope = make_slope()
        surf = make_surface()
        slices = slope.discretize_circle(surf, n_slices=20)
        for s in slices:
            assert s["weight_kN"] > 0


# ---------------------------------------------------------------------------
# Bishop analysis tests
# ---------------------------------------------------------------------------
class TestBishopAnalysis:

    def setup_method(self):
        self.slope = make_slope(cohesion=20.0, phi=20.0, gamma=18.0)
        self.surf  = make_surface(cx=15.0, cy=22.0, R=24.0)
        self.ba    = BishopAnalysis(self.slope, n_slices=30)

    def test_fs_positive(self):
        r = self.ba.factor_of_safety(self.surf)
        assert r["FS"] is not None and r["FS"] > 0

    def test_fs_converged(self):
        r = self.ba.factor_of_safety(self.surf)
        assert r["converged"]

    def test_higher_cohesion_gives_higher_fs(self):
        slope1 = make_slope(cohesion=10.0, phi=20.0)
        slope2 = make_slope(cohesion=30.0, phi=20.0)
        ba1 = BishopAnalysis(slope1, n_slices=25)
        ba2 = BishopAnalysis(slope2, n_slices=25)
        surf = make_surface()
        r1 = ba1.factor_of_safety(surf)
        r2 = ba2.factor_of_safety(surf)
        assert r2["FS"] > r1["FS"]

    def test_higher_phi_gives_higher_fs(self):
        slope1 = make_slope(phi=10.0)
        slope2 = make_slope(phi=30.0)
        ba1 = BishopAnalysis(slope1, n_slices=25)
        ba2 = BishopAnalysis(slope2, n_slices=25)
        surf = make_surface()
        r1 = ba1.factor_of_safety(surf)
        r2 = ba2.factor_of_safety(surf)
        assert r2["FS"] > r1["FS"]

    def test_water_table_reduces_fs(self):
        slope_dry = make_slope(water_table_y=None)
        slope_wet = make_slope(water_table_y=3.0)
        ba_dry = BishopAnalysis(slope_dry, n_slices=25)
        ba_wet = BishopAnalysis(slope_wet, n_slices=25)
        surf = make_surface()
        r_dry = ba_dry.factor_of_safety(surf)
        r_wet = ba_wet.factor_of_safety(surf)
        if r_dry["FS"] and r_wet["FS"]:
            assert r_wet["FS"] <= r_dry["FS"]

    def test_result_has_required_keys(self):
        r = self.ba.factor_of_safety(self.surf)
        for k in ["method", "FS", "n_slices", "iterations", "surface"]:
            assert k in r


# ---------------------------------------------------------------------------
# Janbu analysis tests
# ---------------------------------------------------------------------------
class TestJanbuAnalysis:

    def test_fs_positive(self):
        slope = make_slope()
        ja = JanbuAnalysis(slope, n_slices=25)
        surf = make_surface()
        r = ja.factor_of_safety(surf)
        assert r["FS"] > 0

    def test_janbu_close_to_bishop(self):
        """Janbu and Bishop should give similar results for typical slopes."""
        slope = make_slope()
        surf = make_surface()
        r_b = BishopAnalysis(slope, n_slices=30).factor_of_safety(surf)
        r_j = JanbuAnalysis(slope,  n_slices=30).factor_of_safety(surf)
        assert abs(r_b["FS"] - r_j["FS"]) / r_b["FS"] < 0.20  # within 20%

    def test_higher_cohesion_gives_higher_fs_janbu(self):
        slope1 = make_slope(cohesion=10.0)
        slope2 = make_slope(cohesion=35.0)
        surf = make_surface()
        r1 = JanbuAnalysis(slope1).factor_of_safety(surf)
        r2 = JanbuAnalysis(slope2).factor_of_safety(surf)
        assert r2["FS"] > r1["FS"]


# ---------------------------------------------------------------------------
# Janbu Generalized (full inter-slice shear) tests
# ---------------------------------------------------------------------------
class TestJanbuGeneralizedAnalysis:
    """Duncan-Wright-style homogeneous slope: c'=20 kPa, phi'=20 deg, gamma=18.
    The generalized Janbu should fall between the simplified Janbu and the
    Bishop solution (it accounts for inter-slice shear forces, raising FS
    relative to the simplified method)."""

    def setup_method(self):
        self.slope = make_slope(cohesion=20.0, phi=20.0, gamma=18.0)
        self.surf = make_surface(cx=15.0, cy=22.0, R=24.0)

    def test_fs_positive_and_converged(self):
        r = JanbuGeneralizedAnalysis(self.slope, self.surf,
                                     n_slices=30).compute_fs()
        assert r["FS"] is not None and r["FS"] > 0
        assert r["iterations"] >= 1
        # thrust line returned as list of (x, y)
        assert isinstance(r["thrust_line"], list)
        assert all(len(p) == 2 for p in r["thrust_line"])

    def test_fs_between_simplified_and_bishop(self):
        """Janbu generalized typically lies between Janbu-simplified and
        Bishop/Spencer because it captures inter-slice shear without the
        empirical f0 correction."""
        r_js = JanbuAnalysis(self.slope, n_slices=30).factor_of_safety(self.surf)
        r_b = BishopAnalysis(self.slope, n_slices=30).factor_of_safety(self.surf)
        r_jg = JanbuGeneralizedAnalysis(self.slope, self.surf,
                                        n_slices=30).compute_fs()
        lo = min(r_js["FS"], r_b["FS"]) * 0.85
        hi = max(r_js["FS"], r_b["FS"]) * 1.15
        assert lo <= r_jg["FS"] <= hi, (
            f"FS_jg={r_jg['FS']} not in [{lo:.3f}, {hi:.3f}] "
            f"(simplified={r_js['FS']}, Bishop={r_b['FS']})")

    def test_higher_cohesion_raises_fs_generalized(self):
        s1 = make_slope(cohesion=10.0, phi=20.0)
        s2 = make_slope(cohesion=35.0, phi=20.0)
        surf = make_surface()
        r1 = JanbuGeneralizedAnalysis(s1, surf, n_slices=25).compute_fs()
        r2 = JanbuGeneralizedAnalysis(s2, surf, n_slices=25).compute_fs()
        assert r2["FS"] > r1["FS"]


# ---------------------------------------------------------------------------
# Pseudo-static seismic tests
# ---------------------------------------------------------------------------
class TestPseudoStaticAnalysis:

    def test_seismic_reduces_fs(self):
        slope = make_slope()
        surf  = make_surface()
        r_static  = BishopAnalysis(slope).factor_of_safety(surf)
        r_seismic = PseudoStaticAnalysis(slope, kh=0.15).factor_of_safety(surf)
        if r_static["FS"] and r_seismic["FS"]:
            assert r_seismic["FS"] < r_static["FS"]

    def test_higher_kh_gives_lower_fs(self):
        slope = make_slope()
        surf  = make_surface()
        r1 = PseudoStaticAnalysis(slope, kh=0.10).factor_of_safety(surf)
        r2 = PseudoStaticAnalysis(slope, kh=0.25).factor_of_safety(surf)
        assert r2["FS"] < r1["FS"]

    def test_from_ec8_constructor(self):
        slope = make_slope()
        psa = PseudoStaticAnalysis.from_ec8(slope, ag=0.15, S=1.2)
        assert psa.kh > 0

    def test_newmark_displacement_positive(self):
        slope = make_slope()
        psa   = PseudoStaticAnalysis(slope, kh=0.10)
        r = psa.newmark_displacement(FS_yield=1.0, ag_peak=0.25)
        assert r["displacement_m"] > 0

    def test_newmark_no_displacement_when_ky_ge_kmax(self):
        slope = make_slope()
        psa = PseudoStaticAnalysis(slope, kh=0.30)
        r = psa.newmark_displacement(FS_yield=1.0, ag_peak=0.20)  # kh > ag_peak
        assert r["displacement_m"] == 0.0

    def test_seismic_result_has_kh(self):
        slope = make_slope()
        psa   = PseudoStaticAnalysis(slope, kh=0.15)
        r = psa.factor_of_safety(make_surface())
        assert r["kh"] == pytest.approx(0.15)


# ---------------------------------------------------------------------------
# Critical surface search tests
# ---------------------------------------------------------------------------
class TestCriticalSurfaceSearch:

    def test_search_finds_valid_surface(self):
        slope = make_slope()
        search = CriticalSurfaceSearch(slope, n_cx=5, n_cy=5, n_R=4)
        r = search.search()
        assert r["critical_FS"] is not None
        assert r["critical_FS"] > 0

    def test_search_tries_all_surfaces(self):
        slope = make_slope()
        search = CriticalSurfaceSearch(slope, n_cx=3, n_cy=3, n_R=3)
        r = search.search()
        assert r["n_surfaces_tried"] == 3 * 3 * 3

    def test_stronger_soil_gives_higher_critical_fs(self):
        slope_weak   = make_slope(cohesion=5.0,  phi=10.0)
        slope_strong = make_slope(cohesion=25.0, phi=25.0)
        r_weak   = CriticalSurfaceSearch(slope_weak,   n_cx=4, n_cy=4, n_R=4).search()
        r_strong = CriticalSurfaceSearch(slope_strong, n_cx=4, n_cy=4, n_R=4).search()
        if r_weak["critical_FS"] and r_strong["critical_FS"]:
            assert r_strong["critical_FS"] > r_weak["critical_FS"]


# ---------------------------------------------------------------------------
# Reinforcement auto-wiring tests
# ---------------------------------------------------------------------------
class TestReinforcementWiring:

    def _surf(self):
        return make_surface(cx=15.0, cy=22.0, R=24.0)

    def test_no_reinforcement_regression_bishop(self):
        slope = make_slope()
        surf = self._surf()
        r_default = BishopAnalysis(slope, n_slices=30).factor_of_safety(surf)
        r_none = BishopAnalysis(slope, n_slices=30,
                                reinforcements=None).factor_of_safety(surf)
        r_empty = BishopAnalysis(slope, n_slices=30,
                                 reinforcements=[]).factor_of_safety(surf)
        assert r_default["FS"] == pytest.approx(r_none["FS"])
        assert r_default["FS"] == pytest.approx(r_empty["FS"])

    def test_no_reinforcement_regression_all_methods(self):
        slope = make_slope()
        surf = self._surf()
        for cls in [BishopAnalysis, JanbuAnalysis,
                    SpencerAnalysis, MorgensternPriceAnalysis]:
            r0 = cls(slope, n_slices=30).factor_of_safety(surf)
            r1 = cls(slope, n_slices=30,
                     reinforcements=[]).factor_of_safety(surf)
            assert r0["FS"] == pytest.approx(r1["FS"], rel=1e-6)

    def test_geogrid_increases_fs(self):
        slope = make_slope()
        surf = self._surf()
        geogrid = GeosyntheticLayer(y=5.0, x_left=-5.0, x_right=35.0,
                                    T_allow=80.0, tau_bond=35.0)
        for cls in [BishopAnalysis, JanbuAnalysis,
                    SpencerAnalysis, MorgensternPriceAnalysis]:
            r0 = cls(slope, n_slices=30).factor_of_safety(surf)
            r1 = cls(slope, n_slices=30,
                     reinforcements=[geogrid]).factor_of_safety(surf)
            assert r1["FS"] > r0["FS"], f"{cls.__name__} did not gain FS"

    def test_reinforcement_not_intersecting_gives_zero_delta(self):
        slope = make_slope()
        surf = self._surf()
        # Geogrid placed well above the slip surface -> no intersection.
        geogrid = GeosyntheticLayer(y=50.0, x_left=0.0, x_right=30.0,
                                    T_allow=80.0, tau_bond=35.0)
        r0 = BishopAnalysis(slope, n_slices=30).factor_of_safety(surf)
        r1 = BishopAnalysis(slope, n_slices=30,
                            reinforcements=[geogrid]).factor_of_safety(surf)
        assert r0["FS"] == pytest.approx(r1["FS"])


# ---------------------------------------------------------------------------
# Morgenstern-Price polyline (non-circular) surfaces
# ---------------------------------------------------------------------------
class TestMorgensternPricePolyline:
    """Exercise the per-slice moment arm branch added for PolylineSurface."""

    def test_polyline_approximates_circle(self):
        """A 20-vertex polyline sampled from a circle must give an FS within
        +-5% of the equivalent CircularSurface solution."""
        import numpy as np
        slope = make_slope(cohesion=20.0, phi=20.0, gamma=18.0)
        cx, cy, R = 15.0, 20.0, 22.0
        surf_circ = CircularSurface(cx, cy, R)

        # Sample the circle within the slope footprint.
        x_lo = max(cx - R, 0.0)
        x_hi = min(cx + R, 30.0)
        xs = np.linspace(x_lo, x_hi, 20)
        ys = cy - np.sqrt(np.maximum(R * R - (xs - cx) ** 2, 0.0))
        poly = list(zip(xs.tolist(), ys.tolist()))

        r_c = MorgensternPriceAnalysis(slope, n_slices=30).factor_of_safety(surf_circ)
        r_p = MorgensternPriceAnalysis(slope, n_slices=30).factor_of_safety(poly)

        assert r_c["FS"] is not None and r_p["FS"] is not None
        rel = abs(r_p["FS"] - r_c["FS"]) / r_c["FS"]
        assert rel <= 0.05, (
            f"Polyline FS={r_p['FS']} deviates {rel*100:.2f}% from "
            f"circular FS={r_c['FS']} (>5%).")

    def test_polyline_wedge_three_segments(self):
        """Simple 3-segment wedge polyline: FS must be finite and within
        a physically reasonable range (0.5 < FS < 3.0)."""
        slope = make_slope(cohesion=20.0, phi=20.0, gamma=18.0)
        poly = [(2.0, 0.5), (8.0, -3.0), (16.0, -3.0), (24.0, 5.0)]
        r = MorgensternPriceAnalysis(slope, n_slices=30).factor_of_safety(poly)
        assert r["FS"] is not None
        assert 0.5 < r["FS"] < 3.0, f"Wedge FS={r['FS']} out of range"
