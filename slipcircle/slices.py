"""
slices.py — Reusable slice discretization utilities for limit equilibrium methods.

All LEM methods (Bishop, Janbu, Spencer, Morgenstern-Price) operate on a list of
vertical slices. This module provides a single source of truth so the per-slice
quantities are computed consistently across methods.

A slice is represented as a plain dict with these keys:
    x_mid       midpoint x [m]
    x_left      left boundary x [m]
    x_right     right boundary x [m]
    width       b = x_right - x_left [m]
    height      slice height at midpoint [m]
    y_top       ground surface y at midpoint [m]
    y_bot       failure surface y at midpoint [m]
    alpha_rad   inclination of base [rad] (positive upslope)
    weight_kN   total weight per metre run [kN/m]
    cohesion    c' at base [kPa]
    phi_rad     φ' at base [rad]
    u_kPa       pore pressure at base [kPa]
    layer       name of soil layer at base
"""
from __future__ import annotations
import numpy as np
from math import sqrt, atan2, sin, cos
from typing import Iterable, Sequence

from .geometry import Slope, CircularSurface


def slices_from_circle(slope: Slope, surface: CircularSurface,
                       n_slices: int = 30,
                       water_model=None) -> list[dict]:
    """
    Discretise a circular failure surface into ``n_slices`` vertical slices.

    Parameters
    ----------
    slope : Slope
    surface : CircularSurface
    n_slices : int
    water_model : object with ``u(x, y)`` method, optional
        If supplied, overrides ``slope.pore_pressure`` for u computation.

    Returns
    -------
    list of slice-dicts (see module docstring). Empty list if the circle does
    not intersect the slope above ground level.
    """
    cx, cy, R = surface.cx, surface.cy, surface.R
    xs_ground = np.asarray([p[0] for p in slope.profile])

    # Find horizontal extent where circle is below ground.
    x_lo = max(cx - R, float(xs_ground.min()))
    x_hi = min(cx + R, float(xs_ground.max()))
    if x_hi <= x_lo:
        return []

    edges = np.linspace(x_lo, x_hi, n_slices + 1)
    out: list[dict] = []
    for i in range(n_slices):
        xL, xR = float(edges[i]), float(edges[i + 1])
        xm = 0.5 * (xL + xR)
        dx = xR - xL
        disc = R * R - (xm - cx) ** 2
        if disc <= 0:
            continue
        y_bot = cy - sqrt(disc)
        y_top = slope.ground_y(xm)
        h = y_top - y_bot
        if h <= 1e-6:
            continue

        # Base inclination: sin(alpha) = (xm - cx) / R
        alpha = atan2(xm - cx, cy - y_bot)

        soil = slope.layer_at_depth(y_bot)
        gamma = soil.unit_weight
        W = gamma * h * dx + slope.surcharge * dx

        if water_model is not None:
            u = float(water_model.u(xm, y_bot))
        else:
            u = slope.pore_pressure(xm, y_bot)

        # ru override at soil level
        if soil.ru > 0:
            u = max(u, soil.ru * gamma * h)

        out.append({
            "x_mid": xm,
            "x_left": xL,
            "x_right": xR,
            "width": dx,
            "height": h,
            "y_top": y_top,
            "y_bot": y_bot,
            "alpha_rad": alpha,
            "weight_kN": W,
            "cohesion": soil.cohesion,
            "phi_rad": soil.phi_rad,
            "u_kPa": u,
            "layer": soil.name,
        })
    return out


def slices_from_polyline(slope: Slope,
                          polyline: Sequence[tuple[float, float]],
                          n_slices: int = 30,
                          water_model=None) -> list[dict]:
    """
    Discretise a non-circular (polyline) failure surface.

    ``polyline`` is a list of (x, y) points sorted by increasing x defining the
    slip surface from left exit point to right exit point.
    """
    pts = np.asarray(polyline, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2 or len(pts) < 2:
        raise ValueError("polyline must be a list of >=2 (x, y) points")

    x_lo, x_hi = float(pts[0, 0]), float(pts[-1, 0])
    edges = np.linspace(x_lo, x_hi, n_slices + 1)

    xs = pts[:, 0]
    ys = pts[:, 1]

    out: list[dict] = []
    for i in range(n_slices):
        xL, xR = float(edges[i]), float(edges[i + 1])
        xm = 0.5 * (xL + xR)
        dx = xR - xL
        y_bot = float(np.interp(xm, xs, ys))
        y_top = slope.ground_y(xm)
        h = y_top - y_bot
        if h <= 1e-6:
            continue

        # Local slope of failure surface at xm => alpha
        yL = float(np.interp(xL, xs, ys))
        yR = float(np.interp(xR, xs, ys))
        alpha = atan2(yR - yL, dx)  # rises to the right => positive

        soil = slope.layer_at_depth(y_bot)
        gamma = soil.unit_weight
        W = gamma * h * dx + slope.surcharge * dx
        if water_model is not None:
            u = float(water_model.u(xm, y_bot))
        else:
            u = slope.pore_pressure(xm, y_bot)
        if soil.ru > 0:
            u = max(u, soil.ru * gamma * h)

        out.append({
            "x_mid": xm,
            "x_left": xL,
            "x_right": xR,
            "width": dx,
            "height": h,
            "y_top": y_top,
            "y_bot": y_bot,
            "alpha_rad": alpha,
            "weight_kN": W,
            "cohesion": soil.cohesion,
            "phi_rad": soil.phi_rad,
            "u_kPa": u,
            "layer": soil.name,
        })
    return out
