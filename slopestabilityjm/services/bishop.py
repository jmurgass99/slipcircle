"""Servicio Bishop simplificado.

Funcion pura ``analyze`` que toma ``BishopInput`` y devuelve ``BishopResult``.
"""
from __future__ import annotations
from dataclasses import dataclass

from ..bishop import BishopAnalysis
from ..exceptions import GeometryError
from ._common import SlopeSpec, CircleSpec, build_slope, build_circle


@dataclass(frozen=True)
class BishopInput:
    """Parametros de entrada para Bishop simplificado."""
    slope: SlopeSpec
    circle: CircleSpec
    n_slices: int = 30
    tol: float = 1e-6
    max_iter: int = 100


@dataclass(frozen=True)
class BishopResult:
    """Resultado tipado del analisis Bishop."""
    FS: float
    n_slices: int
    iterations: int
    converged: bool
    M_d_kNm: float
    M_reinforcement_kNm: float
    method: str = "Bishop Simplified"


def analyze(inp: BishopInput) -> BishopResult:
    """Calcula el factor de seguridad con Bishop simplificado.

    Lanza ``GeometryError`` si la superficie no produce dovelas validas.
    """
    slope = build_slope(inp.slope)
    surface = build_circle(inp.circle)
    ba = BishopAnalysis(
        slope=slope, n_slices=inp.n_slices,
        tol=inp.tol, max_iter=inp.max_iter,
    )
    r = ba.factor_of_safety(surface)
    if r.get("FS") is None:
        raise GeometryError(r.get("error", "Superficie invalida"))
    return BishopResult(
        FS=float(r["FS"]),
        n_slices=int(r["n_slices"]),
        iterations=int(r["iterations"]),
        converged=bool(r["converged"]),
        M_d_kNm=float(r["M_d_kNm"]),
        M_reinforcement_kNm=float(r.get("M_reinforcement_kNm", 0.0)),
    )
