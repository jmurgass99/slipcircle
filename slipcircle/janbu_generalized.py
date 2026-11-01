"""
janbu_generalized.py — Janbu's Generalized Procedure of Slices (GPS).

Reference
---------
Janbu, N. (1957). Earth pressure and bearing capacity calculations by
generalized procedure of slices. Proc. 4th ICSMFE, London, Vol. 2, 207–212.
Janbu, N. (1973). Slope stability computations. In Embankment-Dam
Engineering, Casagrande Vol., Hirschfeld & Poulos (eds.), Wiley, 47–86.

Difference with the simplified Janbu (``JanbuAnalysis``)
--------------------------------------------------------
The simplified method neglects the difference of inter-slice shear forces
``ΔX = X_{i+1} - X_i`` and applies a correction factor ``f0`` to compensate.
The *generalized* method solves explicitly for the inter-slice shear forces by
imposing moment equilibrium about the base mid-point of each slice. This
provides the *line of thrust* (point of application of the inter-slice normal
force E along the vertical inter-slice boundaries) and the shear-to-normal
ratio X/E along the surface.

Algorithm (per Janbu 1973, Duncan & Wright 2005, §6.5)
------------------------------------------------------
Two nested iterations:

1.  Outer loop on FS using horizontal force equilibrium:

        FS = ( Σ [ (c'·b + (W(1-kv) - u·b - ΔX) · tanφ') · secα / mα ] + H_reinf )
             / ( Σ W(1-kv) · tanα + kh · Σ W )

    with ``mα = cosα · (1 + tanα · tanφ'/FS)``.

2.  Inner pass to update inter-slice shear forces. For each slice the moment
    equation about the base mid-point is

        X_R = (1/b) · [ E_L · (h_L - b·tanα/2) - E_R · (h_R + b·tanα/2)
                        + W·b/2 ] - X_L                                  (*)

    Re-arranged for the *line of thrust* ``h_t = y_E - y_base``:

        h_t(x_R) = [ E_L · h_t(x_L)
                     - (E_R - E_L) · b · tanα / 2
                     + (X_L + X_R) · b / 2
                     - W · b / 2 · ... ]  / E_R

    The classical Janbu update assumes the line of thrust lies at one third of
    the slice height above the base (a robust empirical fallback for narrow
    slices) and then computes X from the inclination of the thrust line:

        X_i ≈ E_i · tan(α_t,i)   with α_t,i = atan( (h_{t,i+1} - h_{t,i}) / b )

    We use a relaxation factor ``omega`` to damp the X update and improve
    convergence.

Returns the FS, the iterations actually performed, the slice count and the
line of thrust ``[(x_R, y_thrust), ...]``.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from math import sin, cos, tan, radians
from typing import List, Optional, Union

from .geometry import Slope, CircularSurface
from .slices import slices_from_circle, slices_from_polyline
from .reinforcement import GeosyntheticLayer, SoilNail


@dataclass
class JanbuGeneralizedAnalysis:
    """
    Janbu Generalized Procedure of Slices (GPS).

    Parameters
    ----------
    slope : Slope
    surface : CircularSurface | list[(x, y)]
        Failure surface. Either a circular surface or a non-circular polyline.
    n_slices : int
        Number of vertical slices used to discretise the surface.
    kh, kv : float
        Pseudo-static seismic coefficients.
    reinforcements : list of GeosyntheticLayer / SoilNail | None
        Reinforcement elements; horizontal components are added to the
        resisting term (same convention as the simplified Janbu).
    max_iter : int
    tol : float
        Convergence tolerance on FS between outer iterations.
    omega : float
        Relaxation factor for the inter-slice shear update (0 < omega <= 1).
        Smaller values increase stability at the price of more iterations.
    """
    slope: Slope
    surface: Optional[Union[CircularSurface, list]] = None
    n_slices: int = 30
    kh: float = 0.0
    kv: float = 0.0
    reinforcements: Optional[List[Union[GeosyntheticLayer, SoilNail]]] = None
    max_iter: int = 100
    tol: float = 1e-4
    omega: float = 0.5

    # ------------------------------------------------------------------
    def _build_slices(self):
        if isinstance(self.surface, CircularSurface):
            return slices_from_circle(self.slope, self.surface, self.n_slices)
        return slices_from_polyline(self.slope, self.surface, self.n_slices)

    # ------------------------------------------------------------------
    def _reinforcement_H(self) -> float:
        """Horizontal resisting force contributed by reinforcements (kN/m)."""
        if not self.reinforcements or self.surface is None:
            return 0.0
        H = 0.0
        surface = self.surface
        is_circ = isinstance(surface, CircularSurface)
        for r in self.reinforcements:
            if isinstance(r, GeosyntheticLayer):
                xi = r.intersect_circle(surface) if is_circ else r.intersect_polyline(surface)
                if xi is None:
                    continue
                sigma_n = max(0.0, (self.slope.ground_y(xi) - r.y) * 18.0)
                T = r.design_force(xi, sigma_n)
                H += T
            elif isinstance(r, SoilNail):
                hit = r.intersect_circle(surface) if is_circ else r.intersect_polyline(surface)
                if hit is None:
                    continue
                _, _, L_b = hit
                T = r.design_force(L_b)
                a = radians(r.inclination_deg)
                H += T * cos(a)
        return H

    # ------------------------------------------------------------------
    def _fs_from_force_eq(self, slices, X_left, X_right, FS_init,
                          H_reinf: float) -> float:
        """
        Iterate Janbu's horizontal-force equilibrium FS equation with the
        current estimate of inter-slice shear differences ``ΔX = X_R - X_L``.
        """
        FS = max(FS_init, 0.3)
        for _ in range(60):
            num = 0.0
            den = 0.0
            FS_prev = FS
            for sl, XL, XR in zip(slices, X_left, X_right):
                alpha = sl["alpha_rad"]
                W = sl["weight_kN"]
                b = sl["width"]
                c = sl["cohesion"]
                phi = sl["phi_rad"]
                u = sl["u_kPa"]
                dX = XR - XL

                Wv = W * (1.0 - self.kv)
                cos_a = max(cos(alpha), 1e-6)
                sec_a = 1.0 / cos_a

                m_alpha = cos_a * (1.0 + tan(alpha) * tan(phi) / FS)
                if abs(m_alpha) < 1e-6:
                    m_alpha = 1e-6 if m_alpha >= 0 else -1e-6

                N_eff = (Wv - u * b - dX)
                shear_resist = (c * b + N_eff * tan(phi)) * sec_a / m_alpha
                num += shear_resist

                den += Wv * tan(alpha) + self.kh * W

            if abs(den) < 1e-10:
                return FS
            FS = (num + H_reinf) / den
            if FS <= 0 or not np.isfinite(FS):
                return max(FS_prev, 0.1)
            if abs(FS - FS_prev) < 1e-6:
                break
        return FS

    # ------------------------------------------------------------------
    def _update_thrust_and_X(self, slices, FS, X_left, X_right,
                             E_left, E_right):
        """
        Update inter-slice normal forces E_i, the line of thrust y_thrust(x)
        and the inter-slice shear forces X_i.

        E_i are obtained by successive integration of the horizontal force
        equation from left to right:

            E_{i+1} = E_i + W·(1-kv)·tanα + kh·W
                          - [(c'·b + (W(1-kv) - u·b - ΔX)·tanφ')/FS] · secα

        The line of thrust h_t,(i+1) is then obtained from moment equilibrium
        about the slice base mid-point. As a robust fallback (Janbu 1973) we
        anchor the thrust at h/3 above the base whenever |E| -> 0.
        """
        n = len(slices)
        E = np.zeros(n + 1)            # inter-slice normal force at boundaries
        h_t = np.zeros(n + 1)          # height of thrust above base at boundary
        new_X = np.zeros(n + 1)

        E[0] = 0.0
        h_t[0] = 0.0

        for i, sl in enumerate(slices):
            alpha = sl["alpha_rad"]
            W = sl["weight_kN"]
            b = sl["width"]
            c = sl["cohesion"]
            phi = sl["phi_rad"]
            u = sl["u_kPa"]
            h_slice = sl["height"]

            Wv = W * (1.0 - self.kv)
            cos_a = max(cos(alpha), 1e-6)
            sec_a = 1.0 / cos_a
            tan_a = tan(alpha)

            dX = X_right[i] - X_left[i]
            m_alpha = cos_a * (1.0 + tan_a * tan(phi) / FS)
            if abs(m_alpha) < 1e-6:
                m_alpha = 1e-6 if m_alpha >= 0 else -1e-6
            N_eff = (Wv - u * b - dX)
            S = (c * b + N_eff * tan(phi)) / FS    # mobilised base shear

            # Horizontal force equilibrium on slice -> E_{i+1}
            # ΣFx = E_R - E_L + kh·W + (S·cosα - N·sinα) - 0 = 0
            # Using the resultant base force, the change in E across a slice is:
            #   dE = (Wv · tanα - S · secα + kh · W) -- after substitution.
            dE = Wv * tan_a + self.kh * W - S * sec_a + (u * b) * tan_a
            E[i + 1] = E[i] + dE

            # Moment about base mid-point (positive CCW). Lever arms:
            #   E_L acts at height h_t[i] on the left boundary (arm b/2 horiz.)
            #   E_R acts at height h_t[i+1] on the right boundary
            #   X_L, X_R act vertically at b/2 horizontal distance
            #   W at centroid (≈ b/2 lever from base mid-point vertically)
            # Solve for h_t[i+1] from ΣM_base_mid = 0 with smooth Tikhonov-style
            # regularization near E -> 0 (toe / crest) toward Janbu (1973) §3.4 h/3.
            E_reg = max(0.05 * W, 1.0)
            if abs(E[i + 1]) > 1e-6:
                M_lhs = (E[i] * h_t[i]
                         + (X_left[i] + X_right[i]) * b / 2.0
                         - E[i] * b / 2.0 * tan_a
                         - E[i + 1] * b / 2.0 * tan_a)
                ht_raw = M_lhs / E[i + 1]
                w_e = (E[i + 1] ** 2) / (E[i + 1] ** 2 + E_reg ** 2)
                h_t[i + 1] = w_e * ht_raw + (1.0 - w_e) * (h_slice / 3.0)
            else:
                h_t[i + 1] = h_slice / 3.0

            # Clamp the thrust line inside the slice for numerical safety.
            h_t[i + 1] = float(np.clip(h_t[i + 1], 0.0, max(h_slice, 1e-3)))

        # Inter-slice shear from inclination of the line of thrust:
        #   X_i ≈ E_i · d(h_t)/dx · (sign convention)
        # We compute the slope of the thrust line between adjacent slices.
        for i in range(1, n):
            b_left = slices[i - 1]["width"]
            b_right = slices[i]["width"]
            # Average slope of thrust line at the i-th vertical boundary
            dx = 0.5 * (b_left + b_right)
            d_ht = h_t[i + 1] - h_t[i - 1]
            if dx > 1e-9:
                tan_at = float(np.clip(d_ht / (2.0 * dx), -1.5, 1.5))
            else:
                tan_at = 0.0
            new_X[i] = E[i] * tan_at

        # Boundary conditions: no inter-slice force at the two ends.
        new_X[0] = 0.0
        new_X[-1] = 0.0
        return E, h_t, new_X

    # ------------------------------------------------------------------
    def compute_fs(self) -> dict:
        slices = self._build_slices()
        if not slices:
            return {"FS": None, "error": "No valid slices",
                    "method": "Janbu Generalized"}

        n = len(slices)
        X_left = np.zeros(n)
        X_right = np.zeros(n)
        H_reinf = self._reinforcement_H()

        FS = 1.5
        E = np.zeros(n + 1)
        h_t = np.zeros(n + 1)

        converged = False
        it_used = 0
        for it in range(self.max_iter):
            it_used = it + 1
            FS_old = FS
            X_old = X_right.copy()

            # Step 1: solve FS from horizontal force equilibrium for current X
            FS = self._fs_from_force_eq(slices, X_left, X_right, FS, H_reinf)

            # Step 2: update X via moment equilibrium / line of thrust
            E, h_t, X_bound = self._update_thrust_and_X(
                slices, FS, X_left, X_right, None, None)

            # Map boundary X (n+1 values) onto slice left/right X (n each)
            new_X_left = X_bound[:-1]
            new_X_right = X_bound[1:]

            # Relaxation to damp oscillations
            X_left = (1.0 - self.omega) * X_left + self.omega * new_X_left
            X_right = (1.0 - self.omega) * X_right + self.omega * new_X_right

            dFS = abs(FS - FS_old)
            dX = float(np.max(np.abs(X_right - X_old))) if n else 0.0
            if dFS < self.tol and dX < 1e-2:
                converged = True
                break

        # Build the line of thrust as absolute coordinates (x, y_thrust)
        thrust_line = []
        # Left exit point (boundary 0)
        thrust_line.append((float(slices[0]["x_left"]),
                            float(slices[0]["y_bot"]) + float(h_t[0])))
        for i, sl in enumerate(slices):
            x_right = float(sl["x_right"])
            y_base_right = float(sl["y_bot"])  # approx: use slice base y
            thrust_line.append((x_right, y_base_right + float(h_t[i + 1])))

        return {
            "method": "Janbu Generalized",
            "FS": round(float(FS), 4),
            "iterations": it_used,
            "n_slices": n,
            "thrust_line": thrust_line,
            "H_reinforcement_kN": round(H_reinf, 2),
            "converged": converged,
        }

    # Keep API compatible with other LEM classes in the library.
    def factor_of_safety(self, surface=None, **kwargs) -> dict:
        """Alias for :meth:`compute_fs`.

        Accepts an optional ``surface`` argument (and ``kh``/``kv``) so the
        signature matches Spencer/MP and the reinforcement-wiring tests, which
        instantiate the class with ``n_slices`` only and then call
        ``factor_of_safety(surface)``.
        """
        if surface is not None:
            self.surface = surface
        for k in ("kh", "kv"):
            if k in kwargs:
                setattr(self, k, kwargs[k])
        return self.compute_fs()
