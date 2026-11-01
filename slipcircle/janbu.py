"""
janbu.py — Janbu Simplified Method for non-circular failure surfaces.

Reference: Janbu, N. (1954). Application of composite slip surface for
stability analysis. Proc. European Conference on Stability of Earth Slopes.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from math import sin, cos, tan, sqrt, radians
from typing import Optional, List, Union
from .geometry import Slope
from .reinforcement import GeosyntheticLayer, SoilNail


@dataclass
class JanbuAnalysis:
    """
    Janbu Simplified Method — satisfies horizontal force equilibrium.

    Applicable to:
    - Non-circular slip surfaces
    - Layered soils
    - Retaining structures

    Parameters
    ----------
    slope : Slope  slope geometry
    n_slices : int  number of slices
    """
    slope: Slope
    n_slices: int = 20
    reinforcements: Optional[List[Union[GeosyntheticLayer, SoilNail]]] = None

    def factor_of_safety(self, surface, correction_f0: float | None = None) -> dict:
        """
        Compute Janbu FS for a given failure surface.

        Parameters
        ----------
        surface : CircularSurface or list of (x,y) points for non-circular
        correction_f0 : float | None
            Janbu correction factor f0 (accounts for interslice shear forces).
            If None, estimated from d/L ratio of the slip surface.
            Typical range 1.05–1.15 for cohesive soils.

        Janbu Simplified Eq.:
          FS = f0 · Σ [c'·b/cos²α + (W - u·b)·tan(φ')/cos²α] / nα(FS)
                    ─────────────────────────────────────────────────────
                            Σ W·tan(α)

          nα = (1 + tan²(α)) × (1/FS + tan(α)·tan(φ'))  -- simplified to:
          nα = cos(α) × [1 + tan(α)·tan(φ')/FS]  / cos(α)
        """
        # For circular surfaces, use same slice discretization
        from .geometry import CircularSurface
        if isinstance(surface, CircularSurface):
            slices = self.slope.discretize_circle(surface, self.n_slices)
        else:
            # Polilínea [(x, y), ...]: misma discretización que Spencer y
            # Morgenstern-Price. (Hasta 2026-10-31 lanzaba NotImplementedError
            # pese a anunciar superficies no circulares.)
            from .slices import slices_from_polyline
            slices = slices_from_polyline(self.slope, surface, self.n_slices)

        if not slices:
            return {"FS": None, "error": "No valid slices"}

        # Correction factor f0 (Janbu, 1973)
        if correction_f0 is None:
            # Estimate from soil type
            has_friction = any(s["phi_rad"] > 0 for s in slices)
            has_cohesion = any(s["cohesion"] > 0 for s in slices)
            if has_friction and has_cohesion:
                correction_f0 = 1.10
            elif has_friction:
                correction_f0 = 1.05
            else:
                correction_f0 = 1.13  # purely cohesive

        # Driving: Σ W·tan(α)
        H_d = sum(s["weight_kN"] * sin(s["alpha_rad"]) /
                  max(cos(s["alpha_rad"]), 1e-6) for s in slices)

        # Reinforcement: extra horizontal resisting force (constant, FS-indep.).
        H_reinf = 0.0
        if self.reinforcements and isinstance(surface, CircularSurface):
            for r in self.reinforcements:
                if isinstance(r, GeosyntheticLayer):
                    xi = r.intersect_circle(surface)
                    if xi is None:
                        continue
                    sigma_n = max(0.0, (self.slope.ground_y(xi) - r.y) * 18.0)
                    T = r.design_force(xi, sigma_n)
                    H_reinf += T  # already horizontal, per metre run
                elif isinstance(r, SoilNail):
                    hit = r.intersect_circle(surface)
                    if hit is None:
                        continue
                    _, _, L_b = hit
                    T = r.design_force(L_b)
                    a = radians(r.inclination_deg)
                    H_reinf += T * cos(a)  # horizontal component

        # Resistance: iterative
        FS = 1.5
        tol = 1e-5
        for _ in range(100):
            FS_old = FS
            num = 0.0
            for s in slices:
                alpha = s["alpha_rad"]
                W = s["weight_kN"]
                b = s["width"]
                c = s["cohesion"]
                phi = s["phi_rad"]
                u = s["u_kPa"]

                cos_a = max(cos(alpha), 1e-6)
                N_eff = W - u * b

                # Janbu nα denominator term
                n_alpha = cos_a * (1 + tan(alpha) * tan(phi) / FS)
                if abs(n_alpha) < 1e-6:
                    n_alpha = 1e-6

                shear = (c * b + N_eff * tan(phi))
                num += shear / n_alpha

            if abs(H_d) < 1e-10:
                break
            FS = correction_f0 * (num + H_reinf) / H_d
            if abs(FS - FS_old) < tol:
                break

        return {
            "method": "Janbu Simplified",
            "FS": round(FS, 4),
            "correction_f0": correction_f0,
            "n_slices": len(slices),
            "H_reinforcement_kN": round(H_reinf, 2),
        }
