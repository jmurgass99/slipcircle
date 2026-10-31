"""
spencer.py — Spencer's Method (1967).

Spencer's method satisfies BOTH force and moment equilibrium of the sliding
mass by assuming that all inter-slice resultant forces are parallel to a
common direction ``theta`` (in radians).  The two unknowns ``(FS, theta)``
are solved iteratively so that both ``Sum F = 0`` (along theta) and
``Sum M = 0`` (about an arbitrary reference point) are satisfied
simultaneously.

For a circular slip surface the moment equation reduces to the classical
``M = R · Sum [Q_i cos(alpha_i - theta)] = 0`` because all base normals
pass through the circle centre.  For a non-circular (polyline) surface this
simplification does not apply and a rigorous moment balance of the full
sliding mass about an arbitrary reference point is computed slice-by-slice,
including the contribution of every base normal/shear and the line of
action of every parallel inter-slice resultant (acting at slice mid-height
— Spencer's original assumption).

The two-equation 1D-in-FS / 1D-in-theta nested iteration is robust and is
identical to the original implementation for circular surfaces.

Reference: Spencer, E. (1967). A Method of Analysis of the Stability of
Embankments Assuming Parallel Inter-Slice Forces. Géotechnique 17(1), 11-26.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from math import sin, cos, tan, sqrt, pi
from typing import Optional, List, Union, Sequence

from scipy.optimize import brentq

from .geometry import Slope, CircularSurface
from .slices import slices_from_circle, slices_from_polyline
from .reinforcement import GeosyntheticLayer, SoilNail, reinforcement_moment


@dataclass
class SpencerAnalysis:
    slope: Slope
    n_slices: int = 30
    tol: float = 1e-4
    max_iter: int = 100
    reinforcements: Optional[List[Union[GeosyntheticLayer, SoilNail]]] = None

    # ------------------------------------------------------------------
    # Per-slice quantities
    # ------------------------------------------------------------------
    def _slice_Q_NS(self, sl: dict, FS: float, theta: float,
                    kh: float, kv: float) -> tuple[float, float, float]:
        """
        Solve the single-slice equilibrium assuming the resultant of the
        inter-slice forces makes a constant angle ``theta`` with horizontal.

        Returns the inter-slice resultant change ``Q_i`` (acting at angle
        theta) and the resulting base normal ``N_i`` and base shear ``S_i``
        consistent with the chosen FS and theta.
        """
        alpha = sl["alpha_rad"]
        W = sl["weight_kN"]
        b = sl["width"]
        c = sl["cohesion"]
        phi = sl["phi_rad"]
        u = sl["u_kPa"]

        Wv = W * (1.0 - kv)
        Wh = kh * W

        # Spencer's per-slice expression for the inter-slice resultant Q_i
        # (acting along the common direction theta).
        num = ((c * b + (Wv - u * b) * tan(phi)) / FS
               - Wv * sin(alpha) - Wh * cos(alpha))
        denom = cos(alpha - theta) * (1.0 + tan(alpha - theta) * tan(phi) / FS)
        if abs(denom) < 1e-9:
            denom = 1e-9 if denom >= 0 else -1e-9
        Q = num / denom

        # Base normal/shear consistent with Q.  Vertical & horizontal force
        # balance on the slice (the inter-slice resultant change has
        # components Q*cos(theta) horizontal and Q*sin(theta) vertical):
        #   N cos(alpha) + S sin(alpha) =  Wv  +  Q sin(theta)
        #   N sin(alpha) - S cos(alpha) =  Wh  -  Q cos(theta)
        rhs_v = Wv + Q * sin(theta)
        rhs_h = Wh - Q * cos(theta)
        N = rhs_v * cos(alpha) + rhs_h * sin(alpha)
        S = rhs_v * sin(alpha) - rhs_h * cos(alpha)
        return Q, N, S

    # ------------------------------------------------------------------
    # Global equilibrium residuals
    # ------------------------------------------------------------------
    def _equilibria_circular(self, slices, FS, theta, cx, cy, kh, kv):
        """
        Force and moment residuals for a circular surface (reference = centre).

        For a circular surface, Spencer (1967, Eq. 9) shows that the moment
        equation about the centre reduces to

            Sum Q_i · cos(alpha_i - theta) = 0

        (the radius R is a common factor and cancels).  This is independent
        of the line of thrust because all parallel inter-slice resultants
        pair-cancel in moment about the centre.
        """
        sumF = 0.0
        sumM = 0.0
        for sl in slices:
            Q, _, _ = self._slice_Q_NS(sl, FS, theta, kh, kv)
            sumF += Q
            sumM += Q * cos(sl["alpha_rad"] - theta)
        return sumF, sumM

    def _equilibria_polyline(self, slices, FS, theta, x_ref, y_ref, kh, kv):
        """
        Rigorous force/moment residuals about an arbitrary point (x_ref, y_ref)
        for a non-circular surface.  All base reactions (N, S) and inter-slice
        resultants Q_i contribute to the moment.
        """
        sumF = 0.0
        sumM = 0.0
        for sl in slices:
            alpha = sl["alpha_rad"]
            xm = sl["x_mid"]
            y_bot = sl["y_bot"]
            y_mid = 0.5 * (sl["y_top"] + y_bot)
            W = sl["weight_kN"]
            Wv = W * (1.0 - kv)
            Wh = kh * W
            b = sl["width"]
            u = sl["u_kPa"]

            Q, N, S = self._slice_Q_NS(sl, FS, theta, kh, kv)
            sumF += Q

            # 1) Weight at slice centroid (xm, y_mid)
            sumM += Wv * (xm - x_ref)             # vertical down: lever = +(xm-xref)
            sumM += Wh * (y_mid - y_ref)          # horizontal (seismic): lever = +(y-yref)

            # 2) Base normal N (perpendicular to base, into mass: dir = (-sin a, cos a))
            #    point of application: base midpoint (xm, y_bot)
            #    M = N_x*(y_b - y_ref) - N_y*(xm - x_ref)
            Nx = -N * sin(alpha)
            Ny =  N * cos(alpha)
            sumM -= (Nx * (y_bot - y_ref) - Ny * (xm - x_ref))

            # 3) Base shear S (along base, opposing motion: dir = (-cos a, -sin a)
            #    for a sliding mass moving downslope to the left in standard
            #    Spencer convention).  Equivalent to: S resists motion.
            Sx =  S * cos(alpha)
            Sy =  S * sin(alpha)
            sumM -= (Sx * (y_bot - y_ref) - Sy * (xm - x_ref))

            # 4) Inter-slice resultant Q (parallel to theta) acting at mid-height
            #    of slice at (xm, y_mid).
            Qx = Q * cos(theta)
            Qy = Q * sin(theta)
            sumM += (Qx * (y_mid - y_ref) - Qy * (xm - x_ref))

        return sumF, sumM

    # ------------------------------------------------------------------
    def factor_of_safety(self, surface, kh: float = 0.0, kv: float = 0.0) -> dict:
        """Solve (FS, theta) so that both global equilibria are satisfied."""
        if isinstance(surface, CircularSurface):
            slices = slices_from_circle(self.slope, surface, self.n_slices)
            cx, cy = surface.cx, surface.cy
            circular = True
        else:
            slices = slices_from_polyline(self.slope, surface, self.n_slices)
            pts = np.asarray(surface, dtype=float)
            # Reference point: chord midpoint of the polyline (arbitrary; the
            # moment equation is invariant w.r.t. the choice of reference IF
            # force equilibrium is satisfied — choosing a point near the mass
            # improves numerical conditioning).
            cx = float(0.5 * (pts[0, 0] + pts[-1, 0]))
            cy = float(0.5 * (pts[0, 1] + pts[-1, 1]))
            circular = False

        if not slices:
            return {"FS": None, "error": "No valid slices"}

        def residuals(FS, theta):
            if circular:
                return self._equilibria_circular(
                    slices, FS, theta, cx, cy, kh, kv)
            return self._equilibria_polyline(
                slices, FS, theta, cx, cy, kh, kv)

        FS = 1.5
        theta = 0.0
        last = (FS, theta)
        it = 0
        for it in range(self.max_iter):
            # (1) given theta, find FS that zeros the force residual
            try:
                FS = brentq(lambda fs: residuals(fs, theta)[0],
                            0.05, 10.0, xtol=self.tol)
            except ValueError:
                # FIXME: no sign change in the FS bracket — fall back to a
                # damped fixed-point step.  This branch is rarely needed for
                # well-posed problems but can be brittle when the surface is
                # nearly tangent to a weak layer.
                FS = max(0.1, FS * 0.95)

            # (2) given FS, find theta that zeros the moment residual
            try:
                theta = brentq(lambda th: residuals(FS, th)[1],
                               -pi / 3, pi / 3, xtol=self.tol)
            except ValueError:
                theta *= 0.95

            if (abs(FS - last[0]) < self.tol * max(abs(FS), 1.0)
                    and abs(theta - last[1]) < self.tol):
                break
            last = (FS, theta)

        sumF, sumM = residuals(FS, theta)

        # ------------------------------------------------------------------
        # Reinforcement: additive resisting-moment convention (matches Bishop)
        # ------------------------------------------------------------------
        M_reinf = 0.0
        if self.reinforcements and isinstance(surface, CircularSurface):
            M_reinf = float(reinforcement_moment(
                self.reinforcements, surface, self.slope))
            M_d = sum(sl["weight_kN"] * sin(sl["alpha_rad"])
                      for sl in slices) * surface.R
            if abs(M_d) > 1e-9:
                FS = FS + M_reinf / M_d

        return {
            "method": "Spencer",
            "FS": round(FS, 4),
            "theta_deg": round(np.degrees(theta), 3),
            "n_slices": len(slices),
            "iterations": it + 1,
            "residual_F": round(sumF, 4),
            "residual_M": round(sumM, 4),
            "M_reinforcement_kNm": round(M_reinf, 2),
            "converged": (abs(sumF) < 1.0 and abs(sumM) < 1.0),
        }
