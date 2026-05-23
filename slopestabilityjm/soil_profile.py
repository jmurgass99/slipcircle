"""
soil_profile.py — Soil layer and water table definitions.
"""
from __future__ import annotations
from dataclasses import dataclass, field


@dataclass
class SoilLayer:
    """
    Single soil stratum with strength and unit weight parameters.

    Parameters
    ----------
    name : str  descriptive name (e.g. 'Fill', 'Clay', 'Sand')
    cohesion : float  effective cohesion c' [kPa] (0 for sand)
    friction_angle : float  effective internal friction angle φ' [deg]
    unit_weight : float  total unit weight γ [kN/m³]
    sat_unit_weight : float  saturated unit weight γsat [kN/m³] (default = γ)
    ru : float  pore-pressure ratio ru = u/(γ·h) [0–1], 0 = drained

    Notes
    -----
    For undrained analysis (φ=0, c=Su) set friction_angle=0 and cohesion=Su.
    """
    name: str
    cohesion: float           # c' [kPa]
    friction_angle: float     # φ' [deg]
    unit_weight: float        # γ [kN/m³]
    sat_unit_weight: float = -1.0  # γsat — if -1 use unit_weight
    ru: float = 0.0           # pore pressure ratio

    def __post_init__(self) -> None:
        if self.sat_unit_weight < 0:
            self.sat_unit_weight = self.unit_weight
        if not (0 <= self.friction_angle < 90):
            raise ValueError(f"friction_angle must be in [0, 90). Got {self.friction_angle}")
        if self.cohesion < 0:
            raise ValueError("cohesion must be ≥ 0")

    @property
    def phi_rad(self) -> float:
        import math
        return math.radians(self.friction_angle)

    @property
    def tan_phi(self) -> float:
        import math
        return math.tan(self.phi_rad)
