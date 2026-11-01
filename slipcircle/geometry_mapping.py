"""Mapeo semántico -> malla de dovelas (``SliceGeometryMapping``):
``stratum_to_slices`` y ``reinforcement_to_slices`` invariantes bajo remallado.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import math


@dataclass
class SliceGeometryMapping:
    """Mapeo entre estratos/refuerzos semánticos y dovelas discretas."""

    stratum_to_slices: Dict[str, List[int]] = field(default_factory=dict)
    reinforcement_to_slices: Dict[int, List[int]] = field(default_factory=dict)
    slices: List[Dict[str, float]] = field(default_factory=list)


def build_slope_slice_mesh(
    height_m: float = 10.0,
    slope_angle_deg: float = 45.0,
    cx_m: float = 4.0,
    cy_m: float = 14.0,
    radius_m: float = 15.0,
    n_slices: int = 20,
    layer_boundary_depth_m: float = 5.0,
    reinforcement_elevations_m: Tuple[float, ...] = (3.0, 6.0),
) -> SliceGeometryMapping:
    """Discretiza la masa deslizante en ``n_slices`` dovelas verticales y mapea
    qué estrato gobierna la base de cada dovela (``stratum_to_slices``) y qué
    geomallas cortan cada dovela (``reinforcement_to_slices``).
    """
    n_s = max(4, int(n_slices))
    cot_b = 1.0 / math.tan(math.radians(max(10.0, min(80.0, slope_angle_deg))))
    x_crest = height_m * cot_b

    x_left = max(0.0, cx_m - math.sqrt(max(radius_m**2 - cy_m**2, 1.0)))
    x_right = cx_m + math.sqrt(max(radius_m**2 - (cy_m - height_m) ** 2, 1.0))
    dx = (x_right - x_left) / float(n_s)

    stratum_to_slices: Dict[str, List[int]] = {"upper_stratum": [], "lower_stratum": []}
    reinforcement_to_slices: Dict[int, List[int]] = {
        idx: [] for idx in range(len(reinforcement_elevations_m))
    }
    slices: List[Dict[str, float]] = []

    for i in range(n_s):
        xm = x_left + (i + 0.5) * dx
        if xm <= 0.0:
            y_top = 0.0
        elif xm <= x_crest:
            y_top = (xm / max(x_crest, 1e-6)) * height_m
        else:
            y_top = height_m

        val_sq = max(radius_m**2 - (xm - cx_m) ** 2, 0.1)
        y_base = cy_m - math.sqrt(val_sq)
        h_slice = max(0.1, y_top - y_base)
        sin_a = max(-0.95, min(0.95, (xm - cx_m) / radius_m))
        alpha_rad = math.asin(sin_a)

        local_depth = y_top - y_base
        stratum = "upper_stratum" if local_depth <= layer_boundary_depth_m else "lower_stratum"
        stratum_to_slices[stratum].append(i)

        for r_idx, y_reinf in enumerate(reinforcement_elevations_m):
            if y_base <= y_reinf <= y_top:
                reinforcement_to_slices[r_idx].append(i)

        slices.append({
            "index": float(i),
            "x_mid_m": round(xm, 4),
            "width_m": round(dx, 4),
            "y_top_m": round(y_top, 4),
            "y_base_m": round(y_base, 4),
            "height_m": round(h_slice, 4),
            "alpha_rad": round(alpha_rad, 5),
        })

    return SliceGeometryMapping(
        stratum_to_slices=stratum_to_slices,
        reinforcement_to_slices=reinforcement_to_slices,
        slices=slices,
    )
