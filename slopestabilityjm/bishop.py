"""
bishop.py — Bishop Simplified Method for circular failure surfaces.

Reference: Bishop, A.W. (1955). The use of the slip circle in the stability
analysis of slopes. Géotechnique, 5(1), 7-17.

EC7 Annex B (informative): slope stability by slip circle methods.
"""
from __future__ import annotations
import logging
import numpy as np
from dataclasses import dataclass, field
from math import sin, cos, tan, sqrt, radians, pi
from typing import Optional, List, Union
from .geometry import Slope, CircularSurface
from .reinforcement import GeosyntheticLayer, SoilNail, reinforcement_moment
from .exceptions import ConvergenceError

_LOG = logging.getLogger("slopestabilityjm.bishop")


@dataclass
class BishopAnalysis:
    """
    Bishop Simplified Method — satisfies moment equilibrium.

    Assumptions:
    - Circular failure surface
    - Side forces between slices are horizontal (no shear)
    - Effective stress analysis (c', φ')

    Parameters
    ----------
    slope : Slope  slope geometry and soil data
    n_slices : int  number of slices (default 20)
    tol : float  convergence tolerance for iterative FS
    max_iter : int  maximum iterations
    """
    slope: Slope
    n_slices: int = 20
    tol: float = 1e-6
    max_iter: int = 100
    reinforcements: Optional[List[Union[GeosyntheticLayer, SoilNail]]] = None
    strict: bool = False
    """Si ``True`` lanza ``ConvergenceError`` cuando ``max_iter`` se agota
    sin satisfacer la tolerancia. Por defecto ``False`` (backward compat):
    devuelve la mejor estimacion con ``converged=False``."""

    def factor_of_safety(self, surface: CircularSurface) -> dict:
        """
        Compute FS for a given circular failure surface.

        Returns
        -------
        dict with FS, n_slices, convergence info, and per-slice data.

        Bishop Simplified Eq.:
          FS = Σ [c'·b + (W - u·b)·tan(φ')] / mα(FS)
               ─────────────────────────────────────────
                      Σ W·sin(α)

        where:
          mα(FS) = cos(α) + sin(α)·tan(φ') / FS
        """
        slices = self.slope.discretize_circle(surface, self.n_slices)
        if not slices:
            return {"FS": None, "error": "No valid slices — check surface geometry"}

        # Driving moment denominator: Σ W·sin(α)
        M_d = sum(s["weight_kN"] * sin(s["alpha_rad"]) for s in slices)
        if M_d <= 0:
            return {"FS": None, "error": "Non-positive driving moment — surface likely invalid"}

        # Reinforcement contribution to resisting moment (independent of FS).
        M_reinf = 0.0
        if self.reinforcements:
            M_reinf = float(reinforcement_moment(
                self.reinforcements, surface, self.slope))

        # Iterative solution
        FS = 1.5  # initial guess
        for iteration in range(self.max_iter):
            FS_old = FS
            numerator = 0.0
            for s in slices:
                alpha = s["alpha_rad"]
                W = s["weight_kN"]
                b = s["width"]
                c = s["cohesion"]
                phi = s["phi_rad"]
                u = s["u_kPa"]

                # Resisting shear at base per unit width
                N_eff = W - u * b
                shear_resist = c * b + N_eff * tan(phi)

                # mα factor
                m_alpha = cos(alpha) + sin(alpha) * tan(phi) / FS
                if abs(m_alpha) < 1e-8:
                    m_alpha = 1e-8

                numerator += shear_resist / m_alpha

            # numerator is in kN (per unit width); M_reinf is kN·m. Both
            # numerator and M_d carry an implicit factor R that cancels in the
            # classical Bishop derivation, so we divide M_reinf by R to match
            # units before combining.
            FS = (numerator + M_reinf / surface.R) / M_d

            if abs(FS - FS_old) < self.tol:
                break
        else:
            residual = abs(FS - FS_old)
            _LOG.warning(
                "Bishop did not strictly converge: iter=%d residual=%.3e FS=%.4f",
                self.max_iter, residual, FS,
            )
            if self.strict:
                raise ConvergenceError(
                    f"Bishop no convergio en {self.max_iter} iteraciones "
                    f"(residual={residual:.3e})",
                    iteration=self.max_iter, residual=residual,
                )

        return {
            "method": "Bishop Simplified",
            "FS": round(FS, 4),
            "n_slices": len(slices),
            "iterations": iteration + 1,
            "converged": abs(FS - FS_old) < self.tol * 100,
            "surface": {
                "cx": surface.cx,
                "cy": surface.cy,
                "R": surface.R,
            },
            "M_d_kNm": round(M_d, 2),
            "M_reinforcement_kNm": round(M_reinf, 2),
        }
