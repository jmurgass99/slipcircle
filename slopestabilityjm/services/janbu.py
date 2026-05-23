"""Servicio Janbu simplificado."""
from __future__ import annotations
from dataclasses import dataclass

from ..janbu import JanbuAnalysis
from ..exceptions import GeometryError
from ._common import SlopeSpec, CircleSpec, build_slope, build_circle


@dataclass(frozen=True)
class JanbuInput:
    slope: SlopeSpec
    circle: CircleSpec
    n_slices: int = 30
    correction_f0: float | None = None


@dataclass(frozen=True)
class JanbuResult:
    FS: float
    n_slices: int
    correction_f0: float
    H_reinforcement_kN: float
    method: str = "Janbu Simplified"


def analyze(inp: JanbuInput) -> JanbuResult:
    """Calcula FS con Janbu simplificado para superficie circular."""
    slope = build_slope(inp.slope)
    surface = build_circle(inp.circle)
    ja = JanbuAnalysis(slope=slope, n_slices=inp.n_slices)
    r = ja.factor_of_safety(surface, correction_f0=inp.correction_f0)
    if r.get("FS") is None:
        raise GeometryError(r.get("error", "Superficie invalida"))
    return JanbuResult(
        FS=float(r["FS"]),
        n_slices=int(r["n_slices"]),
        correction_f0=float(r["correction_f0"]),
        H_reinforcement_kN=float(r.get("H_reinforcement_kN", 0.0)),
    )
