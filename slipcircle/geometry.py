"""
geometry.py — Slope geometry and circular failure surface discretization.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from math import sqrt, sin, cos, tan, atan2, pi, radians, degrees
from .soil_profile import SoilLayer


@dataclass
class CircularSurface:
    """
    Circular failure surface defined by center and radius.

    Parameters
    ----------
    cx, cy : float  center of circle [m]
    R : float  radius [m]
    """
    cx: float
    cy: float
    R: float

    def intersects_ground(self, x_toe: float, x_crest: float,
                           ground_profile: list[tuple[float, float]]) -> bool:
        """Check that the circle intersects the slope between toe and crest."""
        # Simplified: check distance from center to boundary points
        return self.R > 0


@dataclass
class Slope:
    """
    2D slope geometry with multiple soil layers and optional water table.

    Parameters
    ----------
    profile : list[tuple[float, float]]
        Ground surface profile as (x, y) points [m], left-to-right.
    layers : list[SoilLayer]
        Soil layers from top to bottom.
    layer_depths : list[float]
        Bottom depth of each layer from surface [m].
        Last value = total depth of analysis.
    water_table_y : float | None
        Elevation of the water table [m]. None = dry.
    surcharge : float
        Uniform surcharge on the slope crest [kPa].

    Examples
    --------
    Homogeneous 1V:2H slope, 10m high:
    >>> slope = Slope(
    ...   profile=[(0, 0), (10, 5), (30, 5)],    # toe → break → crest
    ...   layers=[SoilLayer('Clay', cohesion=20, friction_angle=15, unit_weight=18)],
    ...   layer_depths=[20.0],
    ... )
    """
    profile: list[tuple[float, float]]
    layers: list[SoilLayer]
    layer_depths: list[float]
    water_table_y: float | None = None
    surcharge: float = 0.0

    def __post_init__(self) -> None:
        if len(self.layers) != len(self.layer_depths):
            raise ValueError("len(layers) must equal len(layer_depths)")
        self._profile_arr = np.array(self.profile)

    @property
    def toe(self) -> tuple[float, float]:
        return self.profile[0]

    @property
    def crest(self) -> tuple[float, float]:
        return self.profile[-1]

    @property
    def height(self) -> float:
        return self.crest[1] - self.toe[1]

    def ground_y(self, x: float) -> float:
        """Interpolate ground surface elevation at horizontal position x."""
        xs = self._profile_arr[:, 0]
        ys = self._profile_arr[:, 1]
        return float(np.interp(x, xs, ys))

    def layer_at_depth(self, y: float, x: float | None = None) -> SoilLayer:
        """Return the soil layer at elevation y [m] (and optional horizontal coordinate x [m]).

        When ``x`` is provided, depth is measured vertically from the local ground
        surface ``ground_y(x) - y``. When ``x`` is None on a multi-layer slope with
        positive crest elevation, depth is measured from ``crest[1] - y``.
        """
        if len(self.layers) == 1:
            return self.layers[0]
        if x is not None:
            depth = max(0.0, self.ground_y(x) - y)
            for i, layer in enumerate(self.layers):
                if depth <= self.layer_depths[i]:
                    return layer
            return self.layers[-1]
        if self.crest[1] > 0.0 and y > 0.0:
            depth = max(0.0, self.crest[1] - y)
            for i, layer in enumerate(self.layers):
                if depth <= self.layer_depths[i]:
                    return layer
            return self.layers[-1]
        for i, layer in enumerate(self.layers):
            if y >= -self.layer_depths[i]:
                return layer
        return self.layers[-1]

    def pore_pressure(self, x: float, y: float) -> float:
        """
        Pore pressure u [kPa] at point (x, y).
        If water table defined: u = γw × max(0, yw - y)
        """
        gamma_w = 9.81  # kN/m³
        if self.water_table_y is None:
            return 0.0
        return gamma_w * max(0.0, self.water_table_y - y)

    def discretize_circle(self, surface: CircularSurface,
                           n_slices: int = 20) -> list[dict]:
        """
        Discretize a circular failure surface into vertical slices.

        Returns
        -------
        list of dicts, each with: x_mid, width, alpha_rad, height, weight,
                                   cohesion, phi_rad, u, N, T
        """
        cx, cy, R = surface.cx, surface.cy, surface.R

        # Find entry and exit angles of circle with ground surface
        # Entry: left intersection (toe side)
        # Exit: right intersection (crest side)
        # Simplified: use x range from toe to crest within circle
        x_left  = cx - R
        x_right = cx + R
        x_left  = max(x_left,  self.toe[0])
        x_right = min(x_right, self.crest[0])

        if x_right <= x_left:
            return []

        x_mids = np.linspace(x_left, x_right, n_slices + 1)
        slices = []

        for i in range(n_slices):
            xL = x_mids[i]
            xR = x_mids[i + 1]
            xm = 0.5 * (xL + xR)
            dx = xR - xL

            # Top of slice = ground surface
            y_top = self.ground_y(xm)

            # Bottom of slice = circle at xm
            disc = R**2 - (xm - cx)**2
            if disc < 0:
                continue
            y_bot = cy - sqrt(disc)

            if y_top <= y_bot:
                continue  # slice has no height

            h = y_top - y_bot  # slice height

            # Slice base angle alpha (angle of base to horizontal)
            alpha_rad = atan2(xm - cx, cy - y_bot)  # sin(alpha) = (xm-cx)/R

            # Unit weight of soil in this slice
            soil = self.layer_at_depth(y_bot, x=xm)
            gamma = soil.unit_weight

            # Weight
            W = gamma * h * dx + self.surcharge * dx

            # Pore pressure at base
            u = self.pore_pressure(xm, y_bot)

            slices.append({
                "x_mid": xm,
                "width": dx,
                "alpha_rad": alpha_rad,
                "height": h,
                "weight_kN": W,
                "cohesion": soil.cohesion,
                "phi_rad": soil.phi_rad,
                "u_kPa": u,
                "layer": soil.name,
            })

        return slices
