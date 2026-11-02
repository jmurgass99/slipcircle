"""
morgenstern_price.py — Morgenstern-Price (1965) method.

Generalises Spencer by letting the inclination of the inter-slice force vary
along the slip surface with a prescribed function ``f(x)``: at every slice
boundary ``x_j`` the shear and normal inter-slice forces satisfy

    X_j = lambda · f(x_j) · E_j

The two unknowns are ``FS`` and ``lambda``; they are solved so that force and
moment equilibrium are both satisfied.

Force equilibrium is imposed slice by slice in the horizontal AND the vertical
direction: starting from E_0 = 0 at the first boundary, each slice gives a 2×2
linear system for its base normal N_i and the force E on its right boundary;
the force residual is the E left over at the last boundary (it must be zero).
This is the general limit equilibrium formulation (Fredlund & Krahn, 1977).
The moment residual is the moment of the external forces (weight, seismic
force, base normal and base shear) about a reference point; once the force
residual vanishes it does not depend on that point, so a circle (about its
centre) and the same circle sampled as a polyline (about the centre of the
circle through its two ends and its middle, see ``_moment_axis``) give the
same FS. Up to 1.3.1 the method used Spencer's resultant Q_i
with a variable inclination θ_i, which only imposes ΣQ·cos θ = 0: with the
half-sine function the circle gave 1.367 and the polyline 1.354.

References: Morgenstern, N.R. & Price, V.E. (1965). The Analysis of the
Stability of General Slip Surfaces. Géotechnique 15(1), 79-93.
Fredlund, D.G. & Krahn, J. (1977). Comparison of slope stability methods of
analysis. Canadian Geotechnical Journal 14(3), 429-439.
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


def _moment_axis(slices, fallback: tuple[float, float]) -> tuple[float, float]:
    """Eje de momentos de una superficie no circular: centro del círculo que pasa por
    los dos extremos de la superficie y por la base de la dovela central.

    Con el equilibrio de fuerzas cumplido el FS no depende del eje, pero el residuo de
    momentos fuera de la solución sí: respecto al centroide de la masa no es monótono
    en FS y la raíz FS_m(λ) se pierde. Con este eje el residuo se comporta como en un
    círculo (y una polilínea muestreada de un círculo usa su centro). Si los tres
    puntos están alineados, o el centro queda por debajo de la superficie, se usa
    ``fallback``.
    """
    x1, y1 = slices[0]["x_left"], slices[0]["y_bot"]
    x3, y3 = slices[-1]["x_right"], slices[-1]["y_bot"]
    mid = slices[len(slices) // 2]
    x2, y2 = mid["x_mid"], mid["y_bot"]
    d = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    chord = hypot(x3 - x1, y3 - y1)
    if abs(d) < 1e-9 * max(chord, 1.0) ** 2:
        return fallback
    s1, s2, s3 = x1 * x1 + y1 * y1, x2 * x2 + y2 * y2, x3 * x3 + y3 * y3
    ux = (s1 * (y2 - y3) + s2 * (y3 - y1) + s3 * (y1 - y2)) / d
    uy = (s1 * (x3 - x2) + s2 * (x1 - x3) + s3 * (x2 - x1)) / d
    if uy <= y2 or hypot(ux - x2, uy - y2) > 20.0 * max(chord, 1.0):
        return fallback
    return float(ux), float(uy)



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
    def _f_boundaries(self, slices) -> tuple[np.ndarray, np.ndarray]:
        """f(x) en el borde izquierdo y en el derecho de cada dovela, con x normalizada
        sobre [x_izquierdo de la primera, x_derecho de la última]."""
        xl = np.array([sl["x_left"] for sl in slices])
        xr = np.array([sl["x_right"] for sl in slices])
        x0, x1 = float(xl.min()), float(xr.max())
        if x1 - x0 < 1e-9:
            return np.ones_like(xl), np.ones_like(xr)
        fl = np.array([self.f((v - x0) / (x1 - x0)) for v in xl])
        fr = np.array([self.f((v - x0) / (x1 - x0)) for v in xr])
        return fl, fr

    # ------------------------------------------------------------------
    def _residuals(self, slices, fs, lam, f_right, x_ref, y_ref, kh, kv):
        """Residuos (E en el último borde, momento global) para (FS, λ).

        Convenio: la masa desliza hacia −x; sobre la dovela i actúan el peso (0, −W_v),
        la acción sísmica (−W_h, 0) en su centroide, la normal N·(−sen α, cos α) y el
        cortante S·(cos α, sen α) en el punto medio de la base, con
        S = (c·l + (N − u·l)·tan φ)/FS, y las fuerzas entre dovelas (E_izq, X_izq) del
        lado izquierdo y −(E_der, X_der) del derecho, con X = λ·f(x)·E en cada borde.
        Equilibrio horizontal y vertical de la dovela → N y E_der (sistema 2×2).
        """
        E_left = 0.0
        X_left = 0.0
        sumM = 0.0
        for sl, fr in zip(slices, f_right):
            alpha = sl["alpha_rad"]
            sa, ca = sin(alpha), cos(alpha)
            b = sl["width"]
            l = b / max(ca, 1e-9)
            tphi = tan(sl["phi_rad"])
            W = sl["weight_kN"]
            Wv = W * (1.0 - kv)
            Wh = kh * W
            a = (sl["cohesion"] * l - sl["u_kPa"] * l * tphi) / fs
            t = tphi / fs
            # x: N·(−sen α + t·cos α) − E_der = −(a·cos α − W_h + E_izq)
            # y: N·(cos α + t·sen α) − λ·f_der·E_der = −(a·sen α − W_v + X_izq)
            a11, a12, b1 = -sa + t * ca, -1.0, -(a * ca - Wh + E_left)
            a21, a22, b2 = ca + t * sa, -lam * fr, -(a * sa - Wv + X_left)
            det = a11 * a22 - a12 * a21
            if abs(det) < 1e-12:
                det = 1e-12 if det >= 0 else -1e-12
            N = (b1 * a22 - a12 * b2) / det
            E_right = (a11 * b2 - a21 * b1) / det
            S = a + t * N
            xm = sl["x_mid"]
            y_bot = sl["y_bot"]
            y_mid = 0.5 * (sl["y_top"] + y_bot)
            sumM += _moment(0.0, -Wv, xm, y_mid, x_ref, y_ref)
            sumM += _moment(-Wh, 0.0, xm, y_mid, x_ref, y_ref)
            sumM += _moment(-N * sa, N * ca, xm, y_bot, x_ref, y_ref)
            sumM += _moment(S * ca, S * sa, xm, y_bot, x_ref, y_ref)
            E_left = E_right
            X_left = lam * fr * E_right
        return E_left, sumM

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

        f_left, f_right = self._f_boundaries(slices)
        mx, my = (cx, cy) if circular else _moment_axis(slices, (cx, cy))

        def residuals(fs, lam):
            return self._residuals(slices, fs, lam, f_right, mx, my, kh, kv)

        # Esquema de Spencer generalizado: FS_f(λ) y FS_m(λ) por separado y λ tal
        # que coinciden (ver _parallel_solver). El sistema 2×2 de cada dovela es
        # singular en FS = −tan(α_i − θ_der)·tan φ_i, con tan θ_der = λ·f(x_der):
        # el bracket de FS queda por encima de esos polos.
        alphas = np.array([sl["alpha_rad"] for sl in slices])
        tanphis = np.array([tan(sl["phi_rad"]) for sl in slices])

        def fs_lower(lam_: float) -> float:
            thetas = np.arctan(lam_ * f_right)
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
