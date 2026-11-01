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
from math import sin, cos, tan, atan, pi, hypot
from typing import Callable, Optional, List, Union

from scipy.optimize import brentq

from .geometry import Slope, CircularSurface
from .slices import slices_from_circle, slices_from_polyline
from .reinforcement import (
    GeosyntheticLayer,
    SoilNail,
    reinforcement_moment,
    reinforcement_moment_polyline,
)
from ._parallel_solver import solve_parallel


def _moment(fx: float, fy: float, x: float, y: float, x_ref: float, y_ref: float) -> float:
    """Momento (antihorario positivo) de la fuerza (fx, fy) aplicada en (x, y) respecto a (x_ref, y_ref)."""
    return (x - x_ref) * fy - (y - y_ref) * fx



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
        # Misma Q_i de Spencer (1967) Ec. 7 con θ_i variable: base l = b·sec α
        # y normal W·cos α (corregido 2026-10-31; antes c·b y W − u·b).
        l = b / max(cos(alpha), 1e-9)
        num = ((c * l + (Wv * cos(alpha) - Wh * sin(alpha) - u * l) * tan(phi)) / fs
               - Wv * sin(alpha) - Wh * cos(alpha))
        denom = cos(alpha - theta) * (1.0 + tan(alpha - theta) * tan(phi) / fs)
        if abs(denom) < 1e-9:
            denom = 1e-9 if denom >= 0 else -1e-9
        Q = num / denom

        rhs_v = Wv + Q * sin(theta)
        # Equilibrio horizontal con la acción sísmica W_h hacia el pie (−x):
        # −N·sin α + S·cos α = W_h + Q·cos θ. Hasta 2026-10-31 W_h entraba con el
        # signo opuesto en N y S (estabilizante) pero como motriz en Q, y con
        # k_h ≠ 0 el método no convergía (FS_f = 0,90, FS_m = 0,69 en el talud de
        # Duncan y Wright con k_h = 0,10; Bishop pseudo-estático: 0,956).
        rhs_h = -Wh - Q * cos(theta)
        N = rhs_v * cos(alpha) + rhs_h * sin(alpha)
        S = rhs_v * sin(alpha) - rhs_h * cos(alpha)
        return Q, N, S

    # ------------------------------------------------------------------
    def _global_residuals(self, slices, fs, thetas, x_ref, y_ref, kh, kv):
        """Residuos de fuerza y de momento del equilibrio GLOBAL de la masa deslizante.

        Momento respecto a (x_ref, y_ref) de las fuerzas EXTERIORES: peso (0, −W_v) y
        acción sísmica (−W_h, 0) en el centroide de la dovela, y reacciones de la base
        N·(−sin α, cos α) y S·(cos α, sin α) en su punto medio. Las fuerzas entre
        dovelas son interiores y se cancelan en el momento global. Hasta 2026-10-31
        la rama de polilínea sumaba además el momento de Q_i y mezclaba convenios de
        signo entre W y N/S, de modo que una polilínea muestreada de un círculo no
        reproducía el FS del círculo. Con la referencia en el centro de un círculo, N
        pasa por él y la expresión se reduce a R·Σ[S_i − W_i·sin α_i] =
        R·Σ Q_i·cos(α_i − θ_i): la ecuación de Spencer (1967), que la rama circular
        usaba ya. Con acción sísmica el momento de W_h también entra, en las dos ramas.
        """
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
            sumF += Q * cos(th)
            sumM += _moment(0.0, -Wv, xm, y_mid, x_ref, y_ref)
            sumM += _moment(-Wh, 0.0, xm, y_mid, x_ref, y_ref)
            sumM += _moment(-N * sin(alpha), N * cos(alpha), xm, y_bot, x_ref, y_ref)
            sumM += _moment(S * cos(alpha), S * sin(alpha), xm, y_bot, x_ref, y_ref)
        return sumF, sumM

    def _equilibria_circular(self, slices, fs, lam, fvals, cx, cy, kh, kv):
        return self._global_residuals(slices, fs, np.arctan(lam * fvals), cx, cy, kh, kv)

    def _equilibria_polyline(self, slices, fs, lam, fvals, x_ref, y_ref, kh, kv):
        return self._global_residuals(slices, fs, np.arctan(lam * fvals), x_ref, y_ref, kh, kv)

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

        # Esquema de Spencer generalizado: FS_f(λ) y FS_m(λ) por separado y λ tal
        # que coinciden (ver _parallel_solver), con el bracket de FS por encima
        # de los polos FS = −tan(α_i − θ_i)·tan φ_i.
        alphas = np.array([sl["alpha_rad"] for sl in slices])
        tanphis = np.array([tan(sl["phi_rad"]) for sl in slices])

        def fs_lower(lam_: float) -> float:
            thetas = np.arctan(lam_ * fvals)
            with np.errstate(all="ignore"):
                poles = -np.tan(alphas - thetas) * tanphis
            poles = poles[np.isfinite(poles)]
            return float(poles.max(initial=0.0))

        grid = np.arange(-0.5, 2.001, 0.1)
        sol = solve_parallel(residuals, grid, fs_lower, fs_hi=10.0, tol=self.tol)
        FS, lam = sol.FS, sol.p
        it = sol.n_scanned
        if not np.isfinite(FS):
            return {"method": "Morgenstern-Price", "FS": None, "converged": False,
                    "error": "sin raíz en FS para ningún λ de la malla",
                    "n_slices": len(slices)}
        sumF, sumM = residuals(FS, lam)

        M_reinf = 0.0
        if self.reinforcements:
            if isinstance(surface, CircularSurface):
                M_reinf = float(reinforcement_moment(
                    self.reinforcements, surface, self.slope))
                M_d = sum(sl["weight_kN"] * sin(sl["alpha_rad"])
                          for sl in slices) * surface.R
            else:
                M_reinf = float(reinforcement_moment_polyline(
                    self.reinforcements, surface, self.slope, cx, cy))
                r_eff = max(float(np.mean([hypot(sl["x_mid"] - cx, sl["y_bot"] - cy) for sl in slices])), 1.0)
                M_d = sum(sl["weight_kN"] * max(sin(sl["alpha_rad"]), 0.05)
                          for sl in slices) * r_eff
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
            "FS_force_eq": sol.FS_force,
            "FS_moment_eq": sol.FS_moment,
            "M_reinforcement_kNm": round(M_reinf, 2),
            "converged": bool(sol.converged),
        }
