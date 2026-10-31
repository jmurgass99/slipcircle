"""Servicio Spencer (force + moment equilibrium)."""
from __future__ import annotations
from dataclasses import dataclass

from ..spencer import SpencerAnalysis
from ..exceptions import GeometryError, ConvergenceError
from ._common import SlopeSpec, CircleSpec, build_slope, build_circle


@dataclass(frozen=True)
class SpencerInput:
    slope: SlopeSpec
    circle: CircleSpec
    n_slices: int = 30
    kh: float = 0.0
    kv: float = 0.0
    tol: float = 1e-4
    max_iter: int = 100
    strict: bool = False


@dataclass(frozen=True)
class SpencerResult:
    FS: float
    theta_deg: float
    n_slices: int
    iterations: int
    residual_F: float
    residual_M: float
    converged: bool
    method: str = "Spencer"


def analyze(inp: SpencerInput) -> SpencerResult:
    """Calcula (FS, theta) con Spencer.

    Si ``strict=True`` y no se satisface el residuo objetivo se lanza
    ``ConvergenceError``.
    """
    slope = build_slope(inp.slope)
    surface = build_circle(inp.circle)
    sa = SpencerAnalysis(
        slope=slope, n_slices=inp.n_slices,
        tol=inp.tol, max_iter=inp.max_iter,
    )
    r = sa.factor_of_safety(surface, kh=inp.kh, kv=inp.kv)
    if r.get("FS") is None:
        raise GeometryError(r.get("error", "Superficie invalida"))
    if inp.strict and not r["converged"]:
        raise ConvergenceError(
            "Spencer no satisfizo equilibrio simultaneo",
            iteration=int(r["iterations"]),
            residual=float(abs(r["residual_F"]) + abs(r["residual_M"])),
        )
    return SpencerResult(
        FS=float(r["FS"]),
        theta_deg=float(r["theta_deg"]),
        n_slices=int(r["n_slices"]),
        iterations=int(r["iterations"]),
        residual_F=float(r["residual_F"]),
        residual_M=float(r["residual_M"]),
        converged=bool(r["converged"]),
    )
