"""
seismic.py — Pseudo-static seismic slope stability analysis.

Reference:
  EC8-5 §4.1.3: Seismic slope stability
  Seed & Martin (1966): Seismic coefficient method
  Newmark (1965): Sliding block analysis for permanent displacement
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from math import sin, cos, tan, sqrt, pi
from .bishop import BishopAnalysis
from .geometry import Slope, CircularSurface


@dataclass
class PseudoStaticAnalysis:
    """
    Pseudo-static seismic slope stability.

    Adds horizontal and vertical pseudo-static forces:
        Fh = kh × W  (horizontal, destabilizing)
        Fv = kv × W  (vertical, may be ±)

    Parameters
    ----------
    slope : Slope
    kh : float  horizontal seismic coefficient [dimensionless]
    kv : float  vertical seismic coefficient (default = 0.5 × kh)
    n_slices : int
    tol : float  convergence tolerance

    EC8-5 §4.1.3.1 recommendations:
        kh = α · S / (2·r)   where:
            α = ag/g (PGA ratio)
            S = soil factor
            r = 1 for displacement < 30 mm (conservative)
            r = 1.5 for displacement < 20 cm (relaxed)
    """
    slope: Slope
    kh: float
    kv: float = 0.0
    n_slices: int = 20
    tol: float = 1e-6
    max_iter: int = 100

    @classmethod
    def from_ec8(cls, slope: Slope, ag: float, S: float = 1.0,
                 r: float = 1.0, n_slices: int = 20) -> "PseudoStaticAnalysis":
        """
        Create from EC8-5 parameters.

        Parameters
        ----------
        ag : float  PGA [g] (e.g. 0.15 for moderate seismicity)
        S  : float  soil factor (EC8 Table 3.2)
        r  : float  reduction factor (1.0 or 1.5)
        """
        kh = ag * S / (2 * r)
        kv = 0.5 * kh  # typical EC8-5 recommendation
        return cls(slope=slope, kh=kh, kv=kv, n_slices=n_slices)

    def factor_of_safety(self, surface: CircularSurface) -> dict:
        """
        Pseudo-static FS using modified Bishop formulation.

        Modified slice weight: W* = W + Fv = (1 ± kv)W
        Modified driving: add horizontal force moment Fh × arm

        Bishop modified (Terzaghi 1950, Seed 1966):
          FS = Σ [c'·b + (W·(1-kv) - u·b)·tan(φ')] / mα(FS)
               ─────────────────────────────────────────────────
               Σ W·sin(α) + kh·Σ(W·cos(α))  ← horizontal adds to driving
        """
        slices = self.slope.discretize_circle(surface, self.n_slices)
        if not slices:
            return {"FS": None, "error": "No valid slices"}

        # Driving moment
        M_d = sum(
            s["weight_kN"] * sin(s["alpha_rad"])
            + self.kh * s["weight_kN"] * cos(s["alpha_rad"])
            for s in slices
        )
        if M_d <= 0:
            return {"FS": None, "error": "Non-positive driving moment"}

        FS = 1.5
        for iteration in range(self.max_iter):
            FS_old = FS
            numerator = 0.0
            for s in slices:
                alpha = s["alpha_rad"]
                W  = s["weight_kN"] * (1 - self.kv)
                b  = s["width"]
                c  = s["cohesion"]
                phi = s["phi_rad"]
                u  = s["u_kPa"]

                N_eff = W - u * b
                shear_resist = c * b + N_eff * tan(phi)

                m_alpha = cos(alpha) + sin(alpha) * tan(phi) / FS
                if abs(m_alpha) < 1e-8:
                    m_alpha = 1e-8
                numerator += shear_resist / m_alpha

            FS = numerator / M_d
            if abs(FS - FS_old) < self.tol:
                break

        return {
            "method": "Pseudo-static (Bishop modified)",
            "FS": round(FS, 4),
            "kh": self.kh,
            "kv": self.kv,
            "n_slices": len(slices),
            "iterations": iteration + 1,
        }

    @staticmethod
    def newmark_integration(time: np.ndarray, accel_g: np.ndarray,
                             ky_g: float, downslope_only: bool = True) -> dict:
        """
        Direct double integration of the Newmark sliding-block equation.

        Solves the rigid-block-on-inclined-plane model: whenever the input
        acceleration (in units of g) exceeds the yield acceleration ``ky_g``,
        the relative velocity is integrated using

            v(t+dt) = v(t) + (a(t) - ky·g) · dt

        and clamped to ``v >= 0``; displacement is the time integral of v.

        Parameters
        ----------
        time : ndarray   time series [s]
        accel_g : ndarray   ground acceleration parallel to slope [g]
        ky_g : float    yield coefficient [g]
        downslope_only : bool
            If True (default), only positive (downslope) sliding is integrated,
            matching the asymmetric, one-sided Newmark model used in
            slope-stability practice (Jibson 2007).

        Returns
        -------
        dict with arrays of velocity and displacement and the peak slip [m].
        """
        time = np.asarray(time, dtype=float)
        a = np.asarray(accel_g, dtype=float) * 9.81  # m/s2
        ay = ky_g * 9.81

        n = len(time)
        v = np.zeros(n)
        d = np.zeros(n)
        for i in range(1, n):
            dt = time[i] - time[i - 1]
            exc = a[i] - ay  # downslope direction positive
            if downslope_only:
                if v[i - 1] > 0 or exc > 0:
                    v[i] = v[i - 1] + exc * dt
                    if v[i] < 0:
                        v[i] = 0.0
                else:
                    v[i] = 0.0
            else:
                v[i] = v[i - 1] + (a[i] - np.sign(v[i - 1] or exc) * ay) * dt
            d[i] = d[i - 1] + 0.5 * (v[i] + v[i - 1]) * dt

        return {
            "method": "Newmark (direct double integration)",
            "velocity_m_s": v,
            "displacement_m": d,
            "peak_displacement_m": float(d.max()),
            "ky_g": ky_g,
        }

    def newmark_displacement(self, FS_yield: float,
                              ag_peak: float,
                              duration: float = 15.0) -> dict:
        """
        Estimate Newmark sliding block displacement [m].

        Simplified Ambraseys-Menu (1988) empirical formula:
          log(D) = 0.90 + log[(1-ky/kmax)^2.53 · (ky/kmax)^(-1.09)]

        Parameters
        ----------
        FS_yield : float  static FS at yield threshold (≈ FS from pseudo-static = 1.0)
        ag_peak  : float  PGA [g]
        duration : float  strong motion duration [s]

        Returns
        -------
        dict with estimated displacement [m] and acceptability.
        """
        ky = self.kh  # yield acceleration
        kmax = ag_peak
        if kmax <= 0 or ky >= kmax:
            return {"displacement_m": 0.0,
                    "note": "ky >= kmax → no sliding"}

        import math
        ratio = ky / kmax
        log_D = 0.90 + math.log10((1 - ratio)**2.53 * ratio**(-1.09))
        D_cm = 10 ** log_D
        D_m  = D_cm / 100

        # EC8-5 §4.1.3: D < 3 cm → r=2 acceptable; < 30 cm → r=1.5
        acceptable = D_m < 0.30  # 30 cm limit

        return {
            "method": "Newmark (Ambraseys-Menu)",
            "ky": ky,
            "kmax": kmax,
            "ratio_ky_kmax": round(ratio, 4),
            "displacement_m": round(D_m, 4),
            "displacement_cm": round(D_cm, 2),
            "acceptable_EC8": acceptable,
            "note": "< 30 cm limit per EC8-5" if acceptable else "⚠ > 30 cm — redesign"
        }
