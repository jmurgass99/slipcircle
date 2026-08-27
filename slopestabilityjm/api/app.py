"""FastAPI exponiendo los servicios de SlopeStabilityJm."""
from __future__ import annotations

try:
    from fastapi import APIRouter, FastAPI, HTTPException
    from fastapi.middleware.cors import CORSMiddleware
except ImportError as exc:  # pragma: no cover
    from ..exceptions import DependencyError
    raise DependencyError(
        "fastapi requerido: pip install slopestabilityjm[api]"
    ) from exc

from .._brand import NAME as _BRAND_NAME
from .. import __version__
from ..services import bishop, janbu, spencer, slip_search, seismic_slope
from ..services._common import SlopeSpec, CircleSpec, SoilLayerSpec
from ..exceptions import (
    SlopeError, GeometryError, GeotechError, CodeError, ConvergenceError,
)


def _slope_from_dict(d: dict) -> SlopeSpec:
    """Convierte payload JSON ``slope`` en ``SlopeSpec``."""
    layers = tuple(SoilLayerSpec(**L) for L in d["layers"])
    profile = tuple(tuple(p) for p in d["profile"])
    return SlopeSpec(
        profile=profile, layers=layers,
        water_table_y=d.get("water_table_y"),
        surcharge=d.get("surcharge", 0.0),
    )


def _circle_from_dict(d: dict) -> CircleSpec:
    return CircleSpec(cx=float(d["cx"]), cy=float(d["cy"]), R=float(d["R"]))


def _wrap(exc: Exception) -> HTTPException:
    """Mapea excepciones tipadas a codigos HTTP."""
    if isinstance(exc, (GeotechError, GeometryError, CodeError)):
        return HTTPException(422, str(exc))
    if isinstance(exc, ConvergenceError):
        return HTTPException(409, f"{exc} [iter={exc.iteration}, residual={exc.residual:.3e}]")
    if isinstance(exc, SlopeError):
        return HTTPException(400, str(exc))
    return HTTPException(500, f"Internal: {exc}")


router = APIRouter()


@router.get("/")
def root() -> dict:
    """Lista endpoints disponibles."""
    return {
        "name": _BRAND_NAME,
        "version": __version__,
        "endpoints": [
            "/api/slope/bishop",
            "/api/slope/janbu",
            "/api/slope/spencer",
            "/api/slope/critical_slip",
            "/api/slope/seismic",
        ],
    }


@router.get("/api/slope/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


@router.post("/api/slope/bishop")
def api_bishop(payload: dict) -> dict:
    """Bishop simplificado para superficie circular dada."""
    try:
        inp = bishop.BishopInput(
            slope=_slope_from_dict(payload["slope"]),
            circle=_circle_from_dict(payload["circle"]),
            n_slices=int(payload.get("n_slices", 30)),
            tol=float(payload.get("tol", 1e-6)),
            max_iter=int(payload.get("max_iter", 100)),
        )
        r = bishop.analyze(inp)
        return r.__dict__
    except Exception as exc:
        raise _wrap(exc)


@router.post("/api/slope/janbu")
def api_janbu(payload: dict) -> dict:
    """Janbu simplificado para superficie circular."""
    try:
        inp = janbu.JanbuInput(
            slope=_slope_from_dict(payload["slope"]),
            circle=_circle_from_dict(payload["circle"]),
            n_slices=int(payload.get("n_slices", 30)),
            correction_f0=payload.get("correction_f0"),
        )
        return janbu.analyze(inp).__dict__
    except Exception as exc:
        raise _wrap(exc)


@router.post("/api/slope/spencer")
def api_spencer(payload: dict) -> dict:
    """Spencer (force + moment equilibrium)."""
    try:
        inp = spencer.SpencerInput(
            slope=_slope_from_dict(payload["slope"]),
            circle=_circle_from_dict(payload["circle"]),
            n_slices=int(payload.get("n_slices", 30)),
            kh=float(payload.get("kh", 0.0)),
            kv=float(payload.get("kv", 0.0)),
            tol=float(payload.get("tol", 1e-4)),
            max_iter=int(payload.get("max_iter", 100)),
            strict=bool(payload.get("strict", False)),
        )
        return spencer.analyze(inp).__dict__
    except Exception as exc:
        raise _wrap(exc)


@router.post("/api/slope/critical_slip")
def api_critical_slip(payload: dict) -> dict:
    """Grid search Bishop para superficie critica."""
    try:
        inp = slip_search.SlipSearchInput(
            slope=_slope_from_dict(payload["slope"]),
            n_slices=int(payload.get("n_slices", 20)),
            n_cx=int(payload.get("n_cx", 8)),
            n_cy=int(payload.get("n_cy", 8)),
            n_R=int(payload.get("n_R", 6)),
            max_depth_below_toe=float(payload.get("max_depth_below_toe", 0.5)),
            require_toe_exit=bool(payload.get("require_toe_exit", True)),
        )
        return slip_search.analyze(inp).__dict__
    except Exception as exc:
        raise _wrap(exc)


@router.post("/api/slope/seismic")
def api_seismic(payload: dict) -> dict:
    """Pseudo-estatico EC8-5 + opcional Newmark."""
    try:
        inp = seismic_slope.SeismicSlopeInput(
            slope=_slope_from_dict(payload["slope"]),
            circle=_circle_from_dict(payload["circle"]),
            kh=float(payload["kh"]),
            kv=float(payload.get("kv", 0.0)),
            n_slices=int(payload.get("n_slices", 30)),
            ag_peak=payload.get("ag_peak"),
        )
        return seismic_slope.analyze(inp).__dict__
    except Exception as exc:
        raise _wrap(exc)


def create_router() -> APIRouter:
    """Devuelve el APIRouter (factory para uniformidad con otros repos)."""
    return router


# ---------------------------------------------------------------------------
# App standalone
# ---------------------------------------------------------------------------
app = FastAPI(title=f"{_BRAND_NAME} API", version=__version__)
app.add_middleware(
    CORSMiddleware, allow_origins=["*"],
    allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
)
app.include_router(router)
