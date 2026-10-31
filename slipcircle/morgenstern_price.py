"""
morgenstern_price.py — Morgenstern-Price (1965) method.

Generalises Spencer by allowing the inclination of the inter-slice resultant
to vary along the slip surface following a prescribed function ``f(x)``:

    tan(theta_i) = lambda · f(x_i)

The two unknowns are ``FS`` and ``lambda``; they are solved so that both
global force and moment equilibrium are satisfied.

For circular surfaces moments are taken about the circle centre (the
classical simplification — base normals contribute zero).  For non-circular
(polyline) surfaces the moment equation includes the full contribution of
every base normal, base shear and (variable-inclination) inter-slice
resultant acting at slice mid-height.  This replaces the previous
"approximate centroid" formulation, which biased FS on wedge / piece-wise
surfaces.

Reference: Morgenstern, N.R. & Price, V.E. (1965). The Analysis of the
Stability of General Slip Surfaces. Géotechnique 15(1), 79-93.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from math import sin, cos, tan, atan, pi
from typing import Callable, Optional, List, Union

from scipy.optimize import brentq

from .geometry import Slope, CircularSurface
from .slices import slices_from_circle, slices_from_polyline
from .reinforcement import GeosyntheticLayer, SoilNail, reinforcement_moment


def half_sine(x_norm: float) -> float:
    """Half-sine f(x) on the normalised slip-surface coordinate [0, 1]."""
    return float(np.sin(np.pi * x_norm))


def constant_f(x_norm: float) -> float:
    return 1.0


def trapezoidal(x_norm: float) -> float:
    if x_norm < 0.25:
        return x_norm / 0.25
    if x_norm > 0.75:
        return (1.0 - x_norm) / 0.25
    return 1.0


@dataclass
class MorgensternPriceAnalysis:
    slope: Slope
    n_slices: int = 30
    f: Callable[[float], float] = field(default=half_sine)
    tol: float = 1e-4
    max_iter: int = 100
    reinforcements: Optional[List[Union[GeosyntheticLayer, SoilNail]]] = None

    # ------------------------------------------------------------------
    def _f_values(self, slices) -> np.ndarray:
        x = np.array([sl["x_mid"] for sl in slices])
        x0, x1 = x.min(), x.max()
        if x1 - x0 < 1e-9:
            return np.ones_like(x)
        xn = (x - x0) / (x1 - x0)
        return np.array([self.f(v) for v in xn])

    # ------------------------------------------------------------------
    def _slice_Q_NS(self, sl: dict, fs: float, theta: float,
                    kh: float, kv: float) -> tuple[float, float, float]:
        """Per-slice inter-slice resultant Q (at angle theta) and N, S."""
        alpha = sl["alpha_rad"]
        W = sl["weight_kN"]
        b = sl["width"]
        c = sl["cohesion"]
        phi = sl["phi_rad"]
        u = sl["u_kPa"]
        Wv = W * (1.0 - kv)
        Wh = kh * W
        N_eff = Wv - u * b
        num = (c * b + N_eff * tan(phi)) / fs - Wv * sin(alpha) - Wh * cos(alpha)
        denom = cos(alpha - theta) * (1.0 + tan(alpha - theta) * tan(phi) / fs)
        if abs(denom) < 1e-9:
            denom = 1e-9 if denom >= 0 else -1e-9
        Q = num / denom

        rhs_v = Wv + Q * sin(theta)
        rhs_h = Wh - Q * cos(theta)
        N = rhs_v * cos(alpha) + rhs_h * sin(alpha)
        S = rhs_v * sin(alpha) - rhs_h * cos(alpha)
        return Q, N, S

    # ------------------------------------------------------------------
    def _equilibria_circular(self, slices, fs, lam, fvals, cx, cy, kh, kv):
        thetas = np.arctan(lam * fvals)
        sumF = 0.0
        sumM = 0.0
        for sl, th, fv in zip(slices, thetas, fvals):
            Q, _, _ = self._slice_Q_NS(sl, fs, float(th), kh, kv)
            sumF += Q * cos(float(th))
            xm = sl["x_mid"]
            y = 0.5 * (sl["y_top"] + sl["y_bot"])
            sumM += Q * ((xm - cx) * (-sin(float(th)))
                         + (y - cy) * cos(float(th)))
        return sumF, sumM

    def _equilibria_polyline(self, slices, fs, lam, fvals, x_ref, y_ref, kh, kv):
        thetas = np.arctan(lam * fvals)
        sumF = 0.0
        sumM = 0.0
        for sl, th in zip(slices, thetas):
            th = float(th)
            alpha = sl["alpha_rad"]
            xm = sl["x_mid"]
            y_bot = sl["y_bot"]
            y_mid = 0.5 * (sl["y_top"] + y_bot)
            W = sl["weight_kN"]
            Wv = W * (1.0 - kv)
            Wh = kh * W

            Q, N, S = self._slice_Q_NS(sl, fs, th, kh, kv)
            # Horizontal force balance (sum of horizontal projections of Q)
            sumF += Q * cos(th)

            # Moment about (x_ref, y_ref).
            sumM += Wv * (xm - x_ref)
            sumM += Wh * (y_mid - y_ref)

            Nx = -N * sin(alpha); Ny =  N * cos(alpha)
            sumM -= (Nx * (y_bot - y_ref) - Ny * (xm - x_ref))

            Sx =  S * cos(alpha); Sy =  S * sin(alpha)
            sumM -= (Sx * (y_bot - y_ref) - Sy * (xm - x_ref))

            Qx = Q * cos(th); Qy = Q * sin(th)
            sumM += (Qx * (y_mid - y_ref) - Qy * (xm - x_ref))

        return sumF, sumM

    # ------------------------------------------------------------------
    def factor_of_safety(self, surface, kh: float = 0.0, kv: float = 0.0) -> dict:
        if isinstance(surface, CircularSurface):
            slices = slices_from_circle(self.slope, surface, self.n_slices)
            cx, cy = surface.cx, surface.cy
            circular = True
        else:
            slices = slices_from_polyline(self.slope, surface, self.n_slices)
            circular = False
            if slices:
                # Reference point = exact weighted centroid of the sliding
                # mass.  Using W_i*(x_mid, y_mid_of_slice) gives the true
                # centre of gravity, so each slice contributes its own moment
                # arm rather than relying on the approximate midpoint between
                # polyline end-points.
                W_arr  = np.array([sl["weight_kN"] for sl in slices])
                x_arr  = np.array([sl["x_mid"] for sl in slices])
                ym_arr = np.array([0.5 * (sl["y_top"] + sl["y_bot"])
                                   for sl in slices])
                W_sum = float(W_arr.sum())
                if W_sum > 1e-9:
                    cx = float(np.sum(W_arr * x_arr) / W_sum)
                    cy = float(np.sum(W_arr * ym_arr) / W_sum)
                else:
                    pts = np.asarray(surface, dtype=float)
                    cx = float(0.5 * (pts[0, 0] + pts[-1, 0]))
                    cy = float(0.5 * (pts[0, 1] + pts[-1, 1]))
            else:
                pts = np.asarray(surface, dtype=float)
                cx = float(0.5 * (pts[0, 0] + pts[-1, 0]))
                cy = float(0.5 * (pts[0, 1] + pts[-1, 1]))

        if not slices:
            return {"FS": None, "error": "No valid slices"}

        fvals = self._f_values(slices)

        def residuals(fs, lam):
            if circular:
                return self._equilibria_circular(
                    slices, fs, lam, fvals, cx, cy, kh, kv)
            return self._equilibria_polyline(
                slices, fs, lam, fvals, cx, cy, kh, kv)

        FS = 1.5
        lam = 0.0
        last = (FS, lam)
        it = 0
        for it in range(self.max_iter):
            try:
                FS = brentq(lambda fs: residuals(fs, lam)[0],
                            0.05, 10.0, xtol=self.tol)
            except ValueError:
                # FIXME: brentq bracket lost — damped fixed-point fallback.
                FS = max(0.1, FS * 0.9)
            try:
                lam = brentq(lambda lm: residuals(FS, lm)[1],
                             -1.5, 1.5, xtol=self.tol)
            except ValueError:
                lam *= 0.9
            if abs(FS - last[0]) < self.tol and abs(lam - last[1]) < self.tol:
                break
            last = (FS, lam)

        sumF, sumM = residuals(FS, lam)

        M_reinf = 0.0
        if self.reinforcements and isinstance(surface, CircularSurface):
            M_reinf = float(reinforcement_moment(
                self.reinforcements, surface, self.slope))
            M_d = sum(sl["weight_kN"] * sin(sl["alpha_rad"])
                      for sl in slices) * surface.R
            if abs(M_d) > 1e-9:
                FS = FS + M_reinf / M_d

        return {
            "method": "Morgenstern-Price",
            "FS": round(FS, 4),
            "lambda": round(lam, 4),
            "f_function": self.f.__name__,
            "n_slices": len(slices),
            "iterations": it + 1,
            "residual_F": round(sumF, 4),
            "residual_M": round(sumM, 4),
            "M_reinforcement_kNm": round(M_reinf, 2),
            "converged": (abs(sumF) < 1.0 and abs(sumM) < 1.0),
        }
