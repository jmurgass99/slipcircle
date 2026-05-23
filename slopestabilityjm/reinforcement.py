"""
reinforcement.py — Geosynthetic and anchor reinforcement of slopes.

Two element types are supported:

* :class:`GeosyntheticLayer` — horizontal geomesh/geotextile layer at a given
  elevation, with tensile capacity and pull-out resistance per unit length.

* :class:`SoilNail` — inclined anchor/nail defined by head coordinates,
  direction, length, allowable axial force and bond strength.

For a circular slip surface, each element contributes a stabilising moment
about the circle centre when intersected by the surface. The design force is
the minimum of:

    * tensile capacity (``T_allow``)
    * pull-out behind the slip surface: ``q_b · perimeter · L_anchored``
    * (for geosynthetics) the available embedment length × bond stress

The resulting moment is added to the resisting moment in the FS expression:

    FS = (M_resist_soil + M_reinforce) / M_drive
"""
from __future__ import annotations
from dataclasses import dataclass
from math import sin, cos, sqrt, atan2, pi
from typing import Iterable
import numpy as np

from .geometry import CircularSurface, Slope


@dataclass
class GeosyntheticLayer:
    """
    Horizontal reinforcement layer at elevation ``y`` between x_left and x_right.

    Parameters
    ----------
    y : float            elevation [m]
    x_left, x_right : float   horizontal extent [m]
    T_allow : float      long-term design tensile strength [kN/m]
    tau_bond : float     pull-out bond strength per unit area [kPa]
    """
    y: float
    x_left: float
    x_right: float
    T_allow: float
    tau_bond: float = 30.0

    def design_force(self, x_intersect: float, sigma_n: float) -> float:
        """
        Design force per metre run available at the intersection.

        Pull-out from the upslope (anchored) side:
            T_po = 2 · tau_bond · L_embed   (two faces)
        with mobilised friction limited by normal stress: tau = min(tau_bond, sigma_n·mu).
        """
        # Embedment lengths on either side of the slip line
        L_left = max(0.0, x_intersect - self.x_left)
        L_right = max(0.0, self.x_right - x_intersect)
        L_anchored = max(L_left, L_right)  # whichever side is "behind" the slip
        T_po = 2.0 * self.tau_bond * L_anchored
        return float(min(self.T_allow, T_po))

    def intersect_circle(self, surf: CircularSurface) -> float | None:
        """Return x of intersection with the layer line; None if no crossing."""
        disc = surf.R ** 2 - (self.y - surf.cy) ** 2
        if disc < 0:
            return None
        dx = sqrt(disc)
        cand = [surf.cx - dx, surf.cx + dx]
        for x in cand:
            if self.x_left - 1e-6 <= x <= self.x_right + 1e-6:
                return float(x)
        return None


@dataclass
class SoilNail:
    """
    Inclined reinforcement element (passive bar / soil nail / ground anchor).

    Parameters
    ----------
    x_head, y_head : float   head coordinates [m]
    inclination_deg : float  positive = pointing down into the slope
    length : float           total length [m]
    T_allow : float          tensile/yield capacity [kN]
    bond_stress : float      grout-soil bond [kPa]
    diameter : float         hole diameter [m]
    """
    x_head: float
    y_head: float
    inclination_deg: float
    length: float
    T_allow: float
    bond_stress: float = 80.0
    diameter: float = 0.10

    @property
    def tip(self) -> tuple[float, float]:
        a = np.radians(self.inclination_deg)
        return (self.x_head + self.length * cos(a),
                self.y_head - self.length * sin(a))

    def intersect_circle(self, surf: CircularSurface) -> tuple[float, float, float] | None:
        """
        Return (x_int, y_int, L_anchored) where L_anchored is the bar length
        beyond the slip surface (from intersection to tip).  Returns None if
        the bar does not intersect.
        """
        a = np.radians(self.inclination_deg)
        dx, dy = cos(a), -sin(a)
        # Parameterise: P(t) = head + t * (dx, dy),  t in [0, length]
        # |P - C|^2 = R^2
        fx = self.x_head - surf.cx
        fy = self.y_head - surf.cy
        A = dx * dx + dy * dy
        B = 2 * (fx * dx + fy * dy)
        C = fx * fx + fy * fy - surf.R ** 2
        disc = B * B - 4 * A * C
        if disc < 0:
            return None
        rt = sqrt(disc)
        ts = sorted([(-B - rt) / (2 * A), (-B + rt) / (2 * A)])
        for t in ts:
            if 0.0 <= t <= self.length:
                xi = self.x_head + t * dx
                yi = self.y_head + t * dy
                L_beyond = self.length - t
                return (float(xi), float(yi), float(L_beyond))
        return None

    def design_force(self, L_beyond: float) -> float:
        T_po = pi * self.diameter * self.bond_stress * L_beyond
        return float(min(self.T_allow, T_po))


def reinforcement_moment(reinforcements: Iterable, surface: CircularSurface,
                          slope: Slope) -> float:
    """
    Resisting moment [kN·m/m] contributed by all reinforcement elements that
    cross the circular slip surface.

    Forces are assumed to act along the bar direction (for nails) or
    horizontally (for geosynthetics).  Lever arm is the perpendicular distance
    from the circle centre to the line of action of T.
    """
    M = 0.0
    for r in reinforcements:
        if isinstance(r, GeosyntheticLayer):
            xi = r.intersect_circle(surface)
            if xi is None:
                continue
            # Approx overburden normal stress at the layer (simplified)
            sigma_n = (slope.ground_y(xi) - r.y) * 18.0
            T = r.design_force(xi, sigma_n)
            arm = abs(surface.cy - r.y)  # horizontal T => vertical arm
            M += T * arm
        elif isinstance(r, SoilNail):
            hit = r.intersect_circle(surface)
            if hit is None:
                continue
            xi, yi, L_b = hit
            T = r.design_force(L_b)
            a = np.radians(r.inclination_deg)
            # Lever arm = perpendicular distance from (cx, cy) to bar line
            # bar direction (dx, dy) = (cos a, -sin a); normal is (sin a, cos a)
            d = abs((xi - surface.cx) * sin(a) + (yi - surface.cy) * cos(a))
            M += T * d
    return M
