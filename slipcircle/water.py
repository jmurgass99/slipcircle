"""
water.py — Pore pressure models.

Three models are provided:

* ``HorizontalWaterTable``  — flat phreatic line at constant elevation.
* ``PiezometricLine``       — piecewise-linear piezometric surface y_w(x).
* ``RuModel``               — constant ru = u / (γ·h) applied to every slice.

Each model exposes a ``u(x, y)`` method returning the pore pressure [kPa] at
the base of a slice.  ``RuModel`` ignores (x, y) and is intended to be passed
along with the slice height; use ``RuModel.u_from_slice(slice_dict)`` in that
case.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import numpy as np

GAMMA_W = 9.81  # kN/m³


@dataclass
class HorizontalWaterTable:
    """Flat phreatic surface at elevation ``y_w`` [m]."""
    y_w: float
    gamma_w: float = GAMMA_W

    def u(self, x: float, y: float) -> float:
        return self.gamma_w * max(0.0, self.y_w - y)


@dataclass
class PiezometricLine:
    """
    Polyline piezometric surface defined by a list of (x, y_w) points.

    Linear interpolation between vertices; constant extrapolation outside.
    """
    points: Sequence[tuple[float, float]]
    gamma_w: float = GAMMA_W

    def __post_init__(self) -> None:
        arr = np.asarray(self.points, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 2 or len(arr) < 2:
            raise ValueError("points must be a list of >=2 (x, y_w) tuples")
        order = np.argsort(arr[:, 0])
        self._x = arr[order, 0]
        self._y = arr[order, 1]

    def y_w(self, x: float) -> float:
        return float(np.interp(x, self._x, self._y))

    def u(self, x: float, y: float) -> float:
        return self.gamma_w * max(0.0, self.y_w(x) - y)


@dataclass
class RuModel:
    """
    Constant pore-pressure ratio ru = u / (γ · h).

    Use ``u_from_slice`` to compute u from a slice dict because (x, y) alone
    does not give the soil column height.
    """
    ru: float
    gamma: float = 18.0  # used only if u_from_slice is called without slice

    def u(self, x: float, y: float) -> float:
        # Fallback: cannot compute without h; returns 0.
        return 0.0

    def u_from_slice(self, sl: dict) -> float:
        return self.ru * sl["weight_kN"] / max(sl["width"], 1e-6)
