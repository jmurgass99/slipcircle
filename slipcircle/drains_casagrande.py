"""
Diseño de drenes en taludes según Casagrande (idea #12).

Calcula caudal drenado por drenes horizontales (drilled drains) usando
la ecuación de Casagrande para flujo en talud saturado:

    Q_per_drain = k * H * L * factor

donde k es permeabilidad, H altura piezométrica, L longitud dren.
Devuelve número mínimo de drenes para reducir nivel freático a la cota
objetivo en función del espaciamiento.

Pure Python — multi-code agnóstico (físico).
"""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class DrainDesignInput:
    k_m_s: float                # permeabilidad m/s
    slope_height_m: float       # altura del talud
    target_drawdown_m: float    # cuánto se quiere bajar la pluma
    drain_length_m: float       # longitud de cada dren
    drain_spacing_m: float      # espaciamiento horizontal entre drenes
    factor: float = 0.6         # factor empírico Casagrande (0.5-0.8)
    code: str = "EU"


@dataclass(frozen=True)
class DrainDesignResult:
    Q_per_drain_m3_s: float
    n_drains_required: int
    total_Q_m3_day: float
    code: str


def design_drains(inp: DrainDesignInput) -> DrainDesignResult:
    if inp.k_m_s <= 0 or inp.slope_height_m <= 0:
        raise ValueError("k y H > 0 requeridos")
    if inp.target_drawdown_m <= 0:
        raise ValueError("target_drawdown > 0")
    if inp.drain_length_m <= 0 or inp.drain_spacing_m <= 0:
        raise ValueError("drain_length y spacing > 0")
    c = inp.code.upper()
    if c not in ("EU", "EUROCODE", "USA", "FHWA", "NSR", "CO", "AS", "AU", "SG", "SS"):
        raise ValueError(f"Código desconocido: {inp.code}")
    H = inp.slope_height_m
    Q_one = inp.k_m_s * H * inp.drain_length_m * inp.factor  # m3/s
    # Q objetivo: para bajar nivel target_drawdown sobre área talud frontal
    # (modelo simplificado): Q_target = k * drawdown * area_unit
    Q_target_per_m = inp.k_m_s * inp.target_drawdown_m * inp.slope_height_m
    # número de drenes por metro de talud
    n_per_m = Q_target_per_m / max(Q_one, 1e-12)
    # En tramo igual a spacing necesitas ceil(n_per_m * spacing)
    n_drains = max(1, math.ceil(n_per_m * inp.drain_spacing_m))
    total_Q = n_drains * Q_one * 86400.0  # m3/day
    return DrainDesignResult(
        Q_per_drain_m3_s=Q_one,
        n_drains_required=n_drains,
        total_Q_m3_day=total_Q,
        code=inp.code,
    )


def time_to_drawdown_s(k_m_s: float, drawdown_m: float, distance_m: float,
                      porosity: float = 0.3) -> float:
    """Tiempo aprox de descenso piezométrico (Theis simplificado)."""
    if k_m_s <= 0 or distance_m <= 0:
        raise ValueError("k y distance > 0")
    if drawdown_m < 0:
        raise ValueError("drawdown >= 0")
    return (porosity * distance_m * distance_m) / max(k_m_s * max(drawdown_m, 0.01), 1e-12)
