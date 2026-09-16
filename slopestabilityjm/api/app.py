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
    if isinstance(exc, (GeotechError, GeometryError, CodeError, KeyError, TypeError, ValueError)):
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
            "/api/slope/mononobe",
        ],
    }


@router.get("/api/slope/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


@router.post("/api/slope/bishop", summary="Factor de seguridad de taludes mediante método simplificado de Bishop (1955)")
def api_bishop(payload: dict) -> dict:
    """Bishop simplificado para superficie circular dada."""
    from ..services import bishop
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


@router.post("/api/slope/janbu", summary="Estabilidad de taludes por método riguroso de Janbu para superficies no circulares")
def api_janbu(payload: dict) -> dict:
    """Janbu simplificado para superficie circular."""
    from ..services import janbu
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


@router.post("/api/slope/spencer", summary="Factor de seguridad de taludes por método de Spencer (equilibrio completo)")
def api_spencer(payload: dict) -> dict:
    """Spencer (force + moment equilibrium)."""
    from ..services import spencer
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


@router.post("/api/slope/critical_slip", summary="Búsqueda de círculo crítico de deslizamiento mediante optimización")
def api_critical_slip(payload: dict) -> dict:
    """Grid search Bishop para superficie critica."""
    from ..services import slip_search
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


@router.post("/api/slope/seismic", summary="Estabilidad pseudoestática de taludes bajo aceleración sísmica (Eurocódigo 7 / EC8)")
def api_seismic(payload: dict) -> dict:
    """Pseudo-estatico EC8-5 + opcional Newmark."""
    from ..services import seismic_slope
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


@router.post("/api/slope/mononobe", summary="Análisis pseudoestático de cuña sísmica según Mononobe (1926)")
def api_mononobe(payload: dict) -> dict:
    """Evaluación pseudo-estática planar según Mononobe."""
    from ..mononobe import MononobeInput, evaluate_mononobe
    try:
        m_inp = MononobeInput(
            slope_height_m=float(payload["slope_height_m"]),
            slope_angle_deg=float(payload["slope_angle_deg"]),
            gamma_kN_m3=float(payload["gamma_kN_m3"]),
            c_kPa=float(payload["c_kPa"]),
            phi_deg=float(payload["phi_deg"]),
            kh=float(payload.get("kh", 0.15)),
            kv=float(payload.get("kv", 0.0)),
        )
        FS_min = float(payload.get("FS_min", 1.1))
        res = evaluate_mononobe(m_inp, FS_min=FS_min)
        return {
            "FS_static": res.FS_static,
            "FS_seismic": res.FS_seismic,
            "psi_deg": res.psi_deg,
            "is_safe": res.is_safe,
            "FS_min": FS_min,
            "code_reference": "Mononobe (1926) / Eurocode 8 Part 5",
        }
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
