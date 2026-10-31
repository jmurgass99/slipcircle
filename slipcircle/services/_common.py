"""Utilidades comunes a los servicios.

Construye objetos ``Slope`` y ``CircularSurface`` desde estructuras planas
(``SoilLayerSpec``, ``SlopeSpec``, ``CircleSpec``) que son serializables a
JSON sin esfuerzo (REST API, persistencia, snapshot).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence

from ..geometry import Slope, CircularSurface
from ..soil_profile import SoilLayer
from ..exceptions import GeometryError, GeotechError


@dataclass(frozen=True)
class SoilLayerSpec:
    """Especificacion plana de una capa de suelo."""
    name: str
    cohesion: float          # c' [kPa]
    friction_angle: float    # phi' [deg]
    unit_weight: float       # gamma [kN/m3]
    depth_bottom: float      # profundidad inferior desde superficie [m]
    sat_unit_weight: float = -1.0
    ru: float = 0.0


@dataclass(frozen=True)
class SlopeSpec:
    """Especificacion plana de un talud."""
    profile: tuple[tuple[float, float], ...]
    layers: tuple[SoilLayerSpec, ...]
    water_table_y: float | None = None
    surcharge: float = 0.0


@dataclass(frozen=True)
class CircleSpec:
    """Especificacion plana de una superficie circular."""
    cx: float
    cy: float
    R: float


def build_slope(spec: SlopeSpec) -> Slope:
    """Construye ``Slope`` validando la especificacion plana."""
    if len(spec.profile) < 2:
        raise GeometryError("El perfil debe tener al menos 2 puntos (toe, crest).")
    if not spec.layers:
        raise GeotechError("Debe especificarse al menos una capa de suelo.")
    layers = [
        SoilLayer(
            name=L.name, cohesion=L.cohesion,
            friction_angle=L.friction_angle, unit_weight=L.unit_weight,
            sat_unit_weight=L.sat_unit_weight, ru=L.ru,
        )
        for L in spec.layers
    ]
    depths = [L.depth_bottom for L in spec.layers]
    profile = [tuple(p) for p in spec.profile]
    try:
        return Slope(
            profile=profile, layers=layers, layer_depths=depths,
            water_table_y=spec.water_table_y, surcharge=spec.surcharge,
        )
    except ValueError as exc:
        raise GeotechError(str(exc)) from exc


def build_circle(spec: CircleSpec) -> CircularSurface:
    """Construye ``CircularSurface`` con validacion."""
    if spec.R <= 0:
        raise GeometryError(f"Radio debe ser > 0 (R={spec.R}).")
    return CircularSurface(spec.cx, spec.cy, spec.R)
