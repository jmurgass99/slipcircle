"""Servicio de busqueda de superficie critica (grid Bishop)."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any

from ..search import CriticalSurfaceSearch
from ..exceptions import GeometryError
from ._common import SlopeSpec, build_slope


@dataclass(frozen=True)
class SlipSearchInput:
    slope: SlopeSpec
    n_slices: int = 20
    n_cx: int = 8
    n_cy: int = 8
    n_R: int = 6
    max_depth_below_toe: float = 0.5
    require_toe_exit: bool = True


@dataclass(frozen=True)
class SlipSearchResult:
    critical_FS: float
    cx: float
    cy: float
    R: float
    n_surfaces_tried: int
    n_valid: int
    acceptable_EC7: bool
    method: str = "Bishop Simplified (grid search)"


def analyze(inp: SlipSearchInput) -> SlipSearchResult:
    """Encuentra la superficie circular critica por grid search."""
    slope = build_slope(inp.slope)
    search = CriticalSurfaceSearch(
        slope=slope, n_slices=inp.n_slices,
        n_cx=inp.n_cx, n_cy=inp.n_cy, n_R=inp.n_R,
        max_depth_below_toe=inp.max_depth_below_toe,
        require_toe_exit=inp.require_toe_exit,
    )
    r = search.search()
    if r.get("critical_FS") is None:
        raise GeometryError(r.get("error", "No se encontro superficie valida"))
    surf = r["critical_surface"]
    return SlipSearchResult(
        critical_FS=float(r["critical_FS"]),
        cx=float(surf["cx"]),
        cy=float(surf["cy"]),
        R=float(surf["R"]),
        n_surfaces_tried=int(r["n_surfaces_tried"]),
        n_valid=int(r["n_valid"]),
        acceptable_EC7=bool(r["acceptable_EC7"]),
    )
