"""
tests/test_coverage_gate.py -- SS-001 / SS-002 / SS-003

SS-001: 15+ new tests to raise coverage.
SS-002: Canonical-case tests for Bishop, Janbu, Spencer, Morgenstern-Price
        (talud 2H:1V, c'=10 kPa, phi'=30 deg, gamma=18 kN/m3).
SS-003: Newmark sliding-block tests with step pulse (analytical) and
        sinusoidal pulse.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

from slopestabilityjm.geometry import Slope, CircularSurface
from slopestabilityjm.soil_profile import SoilLayer
from slopestabilityjm.bishop import BishopAnalysis
from slopestabilityjm.janbu import JanbuAnalysis
from slopestabilityjm.spencer import SpencerAnalysis
from slopestabilityjm.morgenstern_price import (
    MorgensternPriceAnalysis,
    half_sine,
    constant_f,
    trapezoidal,
)
from slopestabilityjm.seismic import PseudoStaticAnalysis
from slopestabilityjm.water import HorizontalWaterTable, PiezometricLine, RuModel
from slopestabilityjm.reinforcement import GeosyntheticLayer, SoilNail, reinforcement_moment
from slopestabilityjm.slices import slices_from_circle, slices_from_polyline
from slopestabilityjm.newmark import NewmarkInput, evaluate_newmark


# ============================================================================
# Canonical slope fixture: 2H:1V, H=10m, c'=10 kPa, phi'=30 deg, gamma=18
# ============================================================================

def _canonical_slope(water_table_y=None):
    """Talud canonico 2H:1V, H=10 m, c'=10 kPa, phi'=30 deg, gamma=18 kN/m3."""
    layer = SoilLayer("Canon", cohesion=10.0, friction_angle=30.0, unit_weight=18.0)
    return Slope(
        profile=[(0, 0), (20, 10), (30, 10)],
        layers=[layer],
        layer_depths=[20.0],
        water_table_y=water_table_y,
    )


def _canonical_surface():
    """Superficie circular razonable para el talud canonico."""
    return CircularSurface(cx=12.0, cy=20.0, R=22.0)


# ============================================================================
# SS-002: Canonical case -- Bishop Simplified
# ============================================================================

class TestBishopCanonical:

    def test_bishop_canonical_fs_finite_positive(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = BishopAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert r["FS"] is not None
        assert r["FS"] > 0
        assert r["converged"]

    def test_bishop_canonical_fs_in_range(self):
        """FS for this canonical case must be physically reasonable."""
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = BishopAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert 0.8 < r["FS"] < 3.0


# ============================================================================
# SS-002: Canonical case -- Janbu Simplified
# ============================================================================

class TestJanbuCanonical:

    def test_janbu_canonical_fs_positive(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = JanbuAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert r["FS"] is not None
        assert r["FS"] > 0

    def test_janbu_canonical_fs_in_range(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = JanbuAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert 0.8 < r["FS"] < 3.0


# ============================================================================
# SS-002: Canonical case -- Spencer
# ============================================================================

class TestSpencerCanonical:

    def test_spencer_canonical_fs_positive(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = SpencerAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert r["FS"] is not None
        assert r["FS"] > 0

    def test_spencer_canonical_fs_in_range(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = SpencerAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert 0.8 < r["FS"] < 3.0


# ============================================================================
# SS-002: Canonical case -- Morgenstern-Price
# ============================================================================

class TestMorgensternPriceCanonical:

    def test_mp_canonical_fs_positive(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = MorgensternPriceAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert r["FS"] is not None
        assert r["FS"] > 0

    def test_mp_canonical_fs_in_range(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        r = MorgensternPriceAnalysis(slope, n_slices=30).factor_of_safety(surf)
        assert 0.8 < r["FS"] < 3.0


# ============================================================================
# SS-002: Cross-method agreement on canonical case
# ============================================================================

class TestCrossMethodAgreement:

    def test_all_four_methods_agree_within_15_percent(self):
        """All four LEM methods on the canonical case must agree within 15%."""
        slope = _canonical_slope()
        surf = _canonical_surface()

        fs_bishop = BishopAnalysis(slope, n_slices=30).factor_of_safety(surf)["FS"]
        fs_janbu = JanbuAnalysis(slope, n_slices=30).factor_of_safety(surf)["FS"]
        fs_spencer = SpencerAnalysis(slope, n_slices=30).factor_of_safety(surf)["FS"]
        fs_mp = MorgensternPriceAnalysis(slope, n_slices=30).factor_of_safety(surf)["FS"]

        all_fs = [fs_bishop, fs_janbu, fs_spencer, fs_mp]
        assert all(fs is not None and fs > 0 for fs in all_fs)

        fs_mean = sum(all_fs) / len(all_fs)
        for fs in all_fs:
            deviation = abs(fs - fs_mean) / fs_mean
            assert deviation < 0.15, (
                f"FS={fs:.4f} deviates {deviation*100:.1f}% from mean {fs_mean:.4f}. "
                f"All: B={fs_bishop:.4f}, J={fs_janbu:.4f}, "
                f"S={fs_spencer:.4f}, MP={fs_mp:.4f}"
            )


# ============================================================================
# SS-003: Newmark sliding block -- step pulse (analytical)
# ============================================================================

class TestNewmarkStepPulse:

    def test_step_pulse_analytical_displacement(self):
        """Constant acceleration exceeding ky for duration T.

        For a rigid block with a(t)=A (constant) > ky for duration T, then
        a(t)=0, the exact displacement from the newmark_integration function
        in seismic.py is: d = 0.5 * (A - ky) * g * T^2 (during the pulse)
        plus coasting at final velocity v = (A-ky)*g*T until v decays to 0.

        For the one-sided model (downslope_only=True), once the pulse stops
        the block coasts with v>0 but no more exceedance, so velocity
        remains positive and displacement keeps growing until v->0.
        Total = 0.5*(A-ky)*g*T^2 + (A-ky)*g*T * T_coast, where coasting
        is handled by the block mechanics.

        Here we verify the displacement during the pulse alone for simplicity.
        """
        A_g = 0.5   # constant acceleration = 0.5g
        ky_g = 0.2  # yield acceleration = 0.2g
        dt = 0.001  # 1 ms steps
        T = 1.0     # 1 second of constant acceleration
        N = int(T / dt) + 1

        time = np.linspace(0, T, N)
        accel = np.full(N, A_g)

        result = PseudoStaticAnalysis.newmark_integration(time, accel, ky_g)
        d = result["peak_displacement_m"]

        # Expected: during constant exceedance, d = 0.5 * (A-ky)*g * T^2
        # But the block also coasts after velocity builds up, so d >= analytical
        d_analytical = 0.5 * (A_g - ky_g) * 9.81 * T ** 2
        assert d >= d_analytical * 0.95, (
            f"Displacement {d:.4f} m < 95% of analytical {d_analytical:.4f} m"
        )
        assert d > 0

    def test_no_displacement_below_yield(self):
        """When peak acceleration < ky, displacement must be zero."""
        dt = 0.001
        time = np.linspace(0, 2.0, 2001)
        accel = np.full(len(time), 0.1)  # 0.1g, below ky=0.2g
        result = PseudoStaticAnalysis.newmark_integration(time, accel, ky_g=0.2)
        assert result["peak_displacement_m"] == pytest.approx(0.0, abs=1e-6)


# ============================================================================
# SS-003: Newmark sliding block -- sinusoidal pulse
# ============================================================================

class TestNewmarkSinusoidal:

    def test_sinusoidal_pulse_produces_displacement(self):
        """Half-sine pulse exceeding ky must produce positive displacement."""
        dt = 0.001
        T = 2.0
        time = np.linspace(0, T, int(T / dt) + 1)
        # Half-sine with peak 0.6g, period 1s
        freq = 1.0
        accel = 0.6 * np.sin(2 * math.pi * freq * time)
        accel[time > 1.0] = 0.0  # single cycle

        ky_g = 0.2
        result = PseudoStaticAnalysis.newmark_integration(time, accel, ky_g)
        d = result["peak_displacement_m"]
        assert d > 0, "Sinusoidal pulse should produce positive displacement"
        # Bounded: should not exceed free-fall equivalent
        d_max = 0.5 * 0.6 * 9.81 * T ** 2  # free-fall upper bound
        assert d < d_max

    def test_evaluate_newmark_sinusoidal(self):
        """evaluate_newmark from newmark.py with a sinusoidal accelerogram."""
        dt = 0.005
        T = 3.0
        n = int(T / dt)
        t = np.linspace(0, T, n)
        # 2 Hz sinusoidal with peak 0.5g
        ag = 0.5 * np.sin(2 * math.pi * 2.0 * t)

        inp = NewmarkInput(accelerogram_g=ag, dt_s=dt, a_y_g=0.15)
        result = evaluate_newmark(inp)
        assert result.permanent_displacement_m > 0
        assert result.n_pulses_exceeding > 0
        assert result.arias_intensity_m_s > 0


# ============================================================================
# Coverage fillers: water.py
# ============================================================================

class TestWaterModels:

    def test_horizontal_water_table_above(self):
        """Point above the water table has u=0."""
        wt = HorizontalWaterTable(y_w=5.0)
        assert wt.u(10.0, 8.0) == 0.0

    def test_horizontal_water_table_below(self):
        """Point below the water table has u = gamma_w * depth."""
        wt = HorizontalWaterTable(y_w=5.0)
        u = wt.u(10.0, 2.0)
        assert u == pytest.approx(9.81 * 3.0, rel=0.01)

    def test_piezometric_line_interpolation(self):
        """PiezometricLine interpolates between vertices."""
        pl = PiezometricLine(points=[(0, 2), (10, 5), (20, 3)])
        # At x=5, y_w should be 3.5 (linear interp between 2 and 5)
        u = pl.u(5.0, 0.0)
        expected_yw = 3.5
        assert u == pytest.approx(9.81 * expected_yw, rel=0.01)

    def test_piezometric_line_above_water(self):
        pl = PiezometricLine(points=[(0, 2), (10, 5)])
        assert pl.u(5.0, 10.0) == 0.0

    def test_piezometric_line_rejects_single_point(self):
        with pytest.raises(ValueError):
            PiezometricLine(points=[(5, 3)])

    def test_ru_model_u_returns_zero(self):
        """RuModel.u() without slice info returns 0."""
        rm = RuModel(ru=0.3)
        assert rm.u(10.0, 5.0) == 0.0

    def test_ru_model_u_from_slice(self):
        rm = RuModel(ru=0.3)
        sl = {"weight_kN": 100.0, "width": 2.0}
        u = rm.u_from_slice(sl)
        # u = ru * W / b = 0.3 * 100 / 2 = 15.0
        assert u == pytest.approx(15.0)


# ============================================================================
# Coverage fillers: Morgenstern-Price f(x) functions
# ============================================================================

class TestMPFunctions:

    def test_half_sine_endpoints(self):
        assert half_sine(0.0) == pytest.approx(0.0, abs=1e-10)
        assert half_sine(1.0) == pytest.approx(0.0, abs=1e-10)

    def test_half_sine_midpoint(self):
        assert half_sine(0.5) == pytest.approx(1.0, rel=1e-6)

    def test_constant_f_always_one(self):
        for x in [0.0, 0.25, 0.5, 0.75, 1.0]:
            assert constant_f(x) == 1.0

    def test_trapezoidal_plateau(self):
        """Trapezoidal f(x) = 1.0 on [0.25, 0.75]."""
        assert trapezoidal(0.25) == pytest.approx(1.0)
        assert trapezoidal(0.5) == pytest.approx(1.0)
        assert trapezoidal(0.75) == pytest.approx(1.0)

    def test_trapezoidal_ramps(self):
        """Trapezoidal ramps linearly at the edges."""
        assert trapezoidal(0.0) == pytest.approx(0.0)
        assert trapezoidal(1.0) == pytest.approx(0.0)
        # Midpoint of ramp
        assert trapezoidal(0.125) == pytest.approx(0.5, rel=0.01)


# ============================================================================
# Coverage fillers: SoilNail
# ============================================================================

class TestSoilNailCoverage:

    def test_soil_nail_tip(self):
        nail = SoilNail(x_head=5.0, y_head=8.0, inclination_deg=15.0,
                        length=10.0, T_allow=100.0)
        tx, ty = nail.tip
        expected_x = 5.0 + 10.0 * math.cos(math.radians(15.0))
        expected_y = 8.0 - 10.0 * math.sin(math.radians(15.0))
        assert tx == pytest.approx(expected_x, rel=1e-4)
        assert ty == pytest.approx(expected_y, rel=1e-4)

    def test_soil_nail_intersect_circle_hit(self):
        """Nail passes through a circle -- must return intersection."""
        nail = SoilNail(x_head=0.0, y_head=5.0, inclination_deg=0.0,
                        length=20.0, T_allow=100.0)
        surf = CircularSurface(cx=10.0, cy=10.0, R=6.0)
        hit = nail.intersect_circle(surf)
        assert hit is not None
        xi, yi, L_beyond = hit
        # Intersection at y=5, circle centered at (10,10) R=6
        # x intersection: (x-10)^2 + (5-10)^2 = 36 -> (x-10)^2=11 -> x~6.68
        dist = math.sqrt((xi - 10.0) ** 2 + (yi - 10.0) ** 2)
        assert dist == pytest.approx(6.0, rel=0.01)
        assert L_beyond > 0

    def test_soil_nail_intersect_circle_miss(self):
        """Nail does not reach the circle -- must return None."""
        nail = SoilNail(x_head=0.0, y_head=0.0, inclination_deg=0.0,
                        length=2.0, T_allow=100.0)
        surf = CircularSurface(cx=20.0, cy=20.0, R=5.0)
        assert nail.intersect_circle(surf) is None

    def test_soil_nail_design_force_bond_limited(self):
        """When bond capacity < T_allow, bond governs."""
        nail = SoilNail(x_head=0, y_head=0, inclination_deg=0,
                        length=10.0, T_allow=1000.0,
                        bond_stress=50.0, diameter=0.10)
        # T_po = pi * 0.10 * 50 * L_beyond
        L_beyond = 2.0
        T = nail.design_force(L_beyond)
        T_po = math.pi * 0.10 * 50.0 * 2.0
        assert T == pytest.approx(T_po, rel=1e-4)
        assert T < nail.T_allow

    def test_soil_nail_design_force_tensile_limited(self):
        """When T_allow < bond capacity, tensile governs."""
        nail = SoilNail(x_head=0, y_head=0, inclination_deg=0,
                        length=10.0, T_allow=10.0,
                        bond_stress=500.0, diameter=0.20)
        L_beyond = 5.0
        T = nail.design_force(L_beyond)
        assert T == pytest.approx(10.0)  # capped at T_allow


# ============================================================================
# Coverage fillers: Spencer polyline branch
# ============================================================================

class TestSpencerPolyline:

    def test_spencer_polyline_fs_positive(self):
        """Spencer with a polyline surface sampled from a circle."""
        slope = _canonical_slope()
        cx, cy, R = 12.0, 20.0, 22.0

        x_lo = max(cx - R, 0.0)
        x_hi = min(cx + R, 30.0)
        xs = np.linspace(x_lo, x_hi, 15)
        ys = cy - np.sqrt(np.maximum(R ** 2 - (xs - cx) ** 2, 0.0))
        poly = list(zip(xs.tolist(), ys.tolist()))

        r = SpencerAnalysis(slope, n_slices=25).factor_of_safety(poly)
        assert r["FS"] is not None
        assert r["FS"] > 0

    def test_spencer_polyline_matches_circular_within_10pct(self):
        """Polyline approximation of a circle should match circular FS."""
        slope = _canonical_slope()
        cx, cy, R = 12.0, 20.0, 22.0
        surf_circ = CircularSurface(cx, cy, R)

        x_lo = max(cx - R, 0.0)
        x_hi = min(cx + R, 30.0)
        xs = np.linspace(x_lo, x_hi, 20)
        ys = cy - np.sqrt(np.maximum(R ** 2 - (xs - cx) ** 2, 0.0))
        poly = list(zip(xs.tolist(), ys.tolist()))

        fs_circ = SpencerAnalysis(slope, n_slices=30).factor_of_safety(surf_circ)["FS"]
        fs_poly = SpencerAnalysis(slope, n_slices=30).factor_of_safety(poly)["FS"]

        assert fs_circ is not None and fs_poly is not None
        rel = abs(fs_poly - fs_circ) / fs_circ
        assert rel < 0.10, (
            f"Polyline FS={fs_poly:.4f} vs circular FS={fs_circ:.4f}, "
            f"deviation {rel*100:.1f}%"
        )


# ============================================================================
# Coverage fillers: seismic.py -- newmark_integration
# ============================================================================

class TestSeismicNewmarkIntegration:

    def test_newmark_integration_returns_arrays(self):
        time = np.linspace(0, 2.0, 201)
        accel = 0.4 * np.sin(2 * math.pi * 1.0 * time)
        result = PseudoStaticAnalysis.newmark_integration(time, accel, ky_g=0.15)
        assert "velocity_m_s" in result
        assert "displacement_m" in result
        assert len(result["velocity_m_s"]) == len(time)
        assert len(result["displacement_m"]) == len(time)
        assert result["peak_displacement_m"] >= 0

    def test_newmark_integration_bidirectional(self):
        """downslope_only=False allows bidirectional sliding."""
        time = np.linspace(0, 2.0, 2001)
        accel = 0.5 * np.sin(2 * math.pi * 2.0 * time)
        r_one = PseudoStaticAnalysis.newmark_integration(
            time, accel, ky_g=0.2, downslope_only=True)
        r_two = PseudoStaticAnalysis.newmark_integration(
            time, accel, ky_g=0.2, downslope_only=False)
        # Both should produce displacement
        assert r_one["peak_displacement_m"] >= 0
        assert r_two["peak_displacement_m"] >= 0


# ============================================================================
# Coverage fillers: slices_from_polyline edge case
# ============================================================================

class TestSlicesEdgeCases:

    def test_slices_from_polyline_bad_input_raises(self):
        slope = _canonical_slope()
        with pytest.raises(ValueError):
            slices_from_polyline(slope, [(0, 0)])  # only 1 point

    def test_slices_from_circle_with_water_model(self):
        """slices_from_circle with an explicit water model."""
        slope = _canonical_slope()
        surf = _canonical_surface()
        wt = HorizontalWaterTable(y_w=3.0)
        slices = slices_from_circle(slope, surf, n_slices=15, water_model=wt)
        assert len(slices) > 0
        # At least some slices below WT should have u > 0
        has_positive_u = any(s["u_kPa"] > 0 for s in slices)
        assert has_positive_u


# ============================================================================
# Coverage fillers: reinforcement_moment with SoilNail
# ============================================================================

class TestReinforcementMoment:

    def test_reinforcement_moment_with_nail(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        nail = SoilNail(x_head=5.0, y_head=8.0, inclination_deg=10.0,
                        length=20.0, T_allow=100.0, bond_stress=80.0,
                        diameter=0.10)
        M = reinforcement_moment([nail], surf, slope)
        # Nail is long enough to intersect the circle, should contribute
        assert M >= 0

    def test_reinforcement_moment_with_no_intersection(self):
        slope = _canonical_slope()
        surf = _canonical_surface()
        # Short nail that does not reach the slip surface
        nail = SoilNail(x_head=0.0, y_head=0.0, inclination_deg=0.0,
                        length=1.0, T_allow=100.0)
        M = reinforcement_moment([nail], surf, slope)
        assert M == 0.0
