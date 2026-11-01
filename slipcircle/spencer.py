"""
spencer.py — Spencer's Method (1967).

Spencer's method satisfies BOTH force and moment equilibrium of the sliding
mass by assuming that all inter-slice resultant forces are parallel to a
common direction ``theta`` (in radians).  The two unknowns ``(FS, theta)``
are solved iteratively so that both ``Sum F = 0`` (along theta) and
``Sum M = 0`` (about an arbitrary reference point) are satisfied
simultaneously.

The two equations are solved as in Spencer (1967): FS_f(theta) from force
equilibrium (Sum Q_i = 0) and FS_m(theta) from the moment of the EXTERNAL
forces of the sliding mass (weight, seismic force at the slice centroid, base
normal and shear at the base midpoint) about a reference point, with theta the
inclination at which both coincide (see ``_parallel_solver``).  For a circular
surface the reference point is the centre and the moment equation reduces to
``R * Sum [Q_i cos(alpha_i - theta)] = 0``; for a polyline it is the chord
midpoint.  ``converged`` reports whether FS_f == FS_m was found.

Reference: Spencer, E. (1967). A Method of Analysis of the Stability of
Embankments Assuming Parallel Inter-Slice Forces. Géotechnique 17(1), 11-26.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from math import sin, cos, tan, sqrt, pi, hypot
from typing import Optional, List, Union, Sequence

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

        # Spencer (1967) Ec. 7: la resistencia de la base se moviliza sobre su
        # longitud l = b·sec α y la normal del peso es W·cos α (menos la
        # componente sísmica W_h·sin α). Hasta 2026-10-31 el numerador usaba
        # c·b y (W − u·b): el de Bishop, no el de Spencer; en φ = 0 no
        # reproducía la forma cerrada c·R·L_arco/ΣW·(x_m − c_x).
        l = b / max(cos(alpha), 1e-9)
        num = ((c * l + (Wv * cos(alpha) - Wh * sin(alpha) - u * l) * tan(phi)) / FS
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
    # Global equilibrium residuals
    # ------------------------------------------------------------------
    def _global_residuals(self, slices, FS, theta, x_ref, y_ref, kh, kv):
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
        for sl in slices:
            alpha = sl["alpha_rad"]
            xm = sl["x_mid"]
            y_bot = sl["y_bot"]
            y_mid = 0.5 * (sl["y_top"] + y_bot)
            W = sl["weight_kN"]
            Wv = W * (1.0 - kv)
            Wh = kh * W
            Q, N, S = self._slice_Q_NS(sl, FS, theta, kh, kv)
            sumF += Q
            sumM += _moment(0.0, -Wv, xm, y_mid, x_ref, y_ref)
            sumM += _moment(-Wh, 0.0, xm, y_mid, x_ref, y_ref)
            sumM += _moment(-N * sin(alpha), N * cos(alpha), xm, y_bot, x_ref, y_ref)
            sumM += _moment(S * cos(alpha), S * sin(alpha), xm, y_bot, x_ref, y_ref)
        return sumF, sumM

    def _equilibria_circular(self, slices, FS, theta, cx, cy, kh, kv):
        return self._global_residuals(slices, FS, theta, cx, cy, kh, kv)

    def _equilibria_polyline(self, slices, FS, theta, x_ref, y_ref, kh, kv):
        return self._global_residuals(slices, FS, theta, x_ref, y_ref, kh, kv)

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

        # Esquema de Spencer: FS_f(θ) y FS_m(θ) por separado y θ tal que coinciden
        # (ver _parallel_solver). El bracket inferior de FS queda por encima de
        # los polos FS = −tan(α_i − θ)·tan φ_i de cada dovela.
        alphas = np.array([sl["alpha_rad"] for sl in slices])
        tanphis = np.array([tan(sl["phi_rad"]) for sl in slices])

        def fs_lower(theta: float) -> float:
            with np.errstate(all="ignore"):
                poles = -np.tan(alphas - theta) * tanphis
            poles = poles[np.isfinite(poles)]
            return float(poles.max(initial=0.0))

        grid = np.radians(np.arange(-10.0, 60.01, 2.5))
        sol = solve_parallel(residuals, grid, fs_lower, fs_hi=10.0, tol=self.tol)
        FS, theta = sol.FS, sol.p
        it = sol.n_scanned
        if not np.isfinite(FS):
            return {"method": "Spencer", "FS": None, "converged": False,
                    "error": "sin raíz en FS para ningún θ de la malla",
                    "n_slices": len(slices)}
        sumF, sumM = residuals(FS, theta)

        # ------------------------------------------------------------------
        # Reinforcement: additive resisting-moment convention (matches Bishop)
        # ------------------------------------------------------------------
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
            "method": "Spencer",
            "FS": round(FS, 4),
            "theta_deg": round(np.degrees(theta), 3),
            "n_slices": len(slices),
            "iterations": it + 1,
            "residual_F": round(sumF, 4),
            "residual_M": round(sumM, 4),
            "FS_force_eq": sol.FS_force,
            "FS_moment_eq": sol.FS_moment,
            "M_reinforcement_kNm": round(M_reinf, 2),
            "converged": bool(sol.converged),
        }
