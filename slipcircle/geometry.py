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
    layer_reference : str
        ``"surface"`` (default): ``layer_depths`` are measured down from the local
        ground surface, so interfaces run parallel to the profile.
        ``"elevation"``: ``layer_depths`` are depths below the datum y = 0, so the
        bottom of layer i is the horizontal line y = -layer_depths[i] (a negative
        value puts it above the datum). Use it for horizontal strata under a
        sloping surface, e.g. an embankment on soft ground: ``layer_depths=[0, 20]``.

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
    layer_reference: str = "surface"

    def __post_init__(self) -> None:
        if len(self.layers) != len(self.layer_depths):
            raise ValueError("len(layers) must equal len(layer_depths)")
        if self.layer_reference not in ("surface", "elevation"):
            raise ValueError("layer_reference must be 'surface' or 'elevation'")
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
        if self.layer_reference == "elevation":
            for i, layer in enumerate(self.layers):
                if y >= -self.layer_depths[i]:
                    return layer
            return self.layers[-1]
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

    def circle_ground_extent(self, surface: CircularSurface,
                             n_probe: int = 4001) -> tuple[float, float] | None:
        """Horizontal range [x_lo, x_hi] where the circle runs below the ground surface.

        The ground beyond the first and last profile points continues horizontally
        (``ground_y`` clamps), so a deep circle that exits beyond the toe keeps the
        soil wedge on that side instead of being cut at the end of the profile.
        Returns None if the circle never goes below ground.
        """
        cx, cy, R = surface.cx, surface.cy, surface.R
        xs = np.linspace(cx - R, cx + R, n_probe)
        yb = cy - np.sqrt(np.maximum(R * R - (xs - cx) ** 2, 0.0))
        yt = np.interp(xs, self._profile_arr[:, 0], self._profile_arr[:, 1])
        below = np.nonzero(yt - yb > 1e-9)[0]
        if below.size == 0:
            return None
        step = xs[1] - xs[0]
        x_lo = max(cx - R, float(xs[below[0]]) - step)
        x_hi = min(cx + R, float(xs[below[-1]]) + step)
        return x_lo, x_hi

    def surcharge_applies(self, x: float) -> bool:
        """Surcharge acts only within the defined profile, not on its horizontal extensions."""
        return float(self._profile_arr[0, 0]) <= x <= float(self._profile_arr[-1, 0])

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
        # x range where the circle is below ground (ground extended horizontally
        # beyond the profile ends, see circle_ground_extent)
        extent = self.circle_ground_extent(surface)
        if extent is None:
            return []
        x_left, x_right = extent

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
            W = gamma * h * dx + (self.surcharge * dx if self.surcharge_applies(xm) else 0.0)

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
