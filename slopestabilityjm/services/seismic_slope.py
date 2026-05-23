"""Servicio sismico pseudo-estatico para taludes (EC8-5)."""
from __future__ import annotations
from dataclasses import dataclass

from ..seismic import PseudoStaticAnalysis
from ..exceptions import GeometryError, CodeError
from ._common import SlopeSpec, CircleSpec, build_slope, build_circle


@dataclass(frozen=True)
class SeismicSlopeInput:
    slope: SlopeSpec
    circle: CircleSpec
    kh: float
    kv: float = 0.0
    n_slices: int = 30
    ag_peak: float | None = None    # opcional: si se da, calcula desplazamiento Newmark


@dataclass(frozen=True)
class SeismicSlopeResult:
    FS: float
    kh: float
    kv: float
    n_slices: int
    iterations: int
    displacement_m: float | None
    acceptable_EC8: bool | None
    method: str = "Pseudo-static (Bishop modified)"


def analyze(inp: SeismicSlopeInput) -> SeismicSlopeResult:
    """Analisis pseudo-estatico + opcional Newmark."""
    if inp.kh < 0:
        raise CodeError(f"kh debe ser >= 0 (kh={inp.kh}).")
    slope = build_slope(inp.slope)
    surface = build_circle(inp.circle)
    psa = PseudoStaticAnalysis(
        slope=slope, kh=inp.kh, kv=inp.kv, n_slices=inp.n_slices,
    )
    r = psa.factor_of_safety(surface)
    if r.get("FS") is None:
        raise GeometryError(r.get("error", "Superficie invalida"))

    disp = None
    acceptable = None
    if inp.ag_peak is not None:
        nm = psa.newmark_displacement(FS_yield=1.0, ag_peak=inp.ag_peak)
        disp = float(nm["displacement_m"])
        acceptable = bool(nm.get("acceptable_EC8", disp < 0.30))

    return SeismicSlopeResult(
        FS=float(r["FS"]),
        kh=float(r["kh"]),
        kv=float(r["kv"]),
        n_slices=int(r["n_slices"]),
        iterations=int(r["iterations"]),
        displacement_m=disp,
        acceptable_EC8=acceptable,
    )
