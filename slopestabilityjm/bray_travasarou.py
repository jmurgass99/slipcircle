"""
Bray-Travasarou (2007) — modelo predictivo de desplazamiento permanente
co-sísmico en taludes/presas, basado en Newmark con T_s (período del
talud) como variable adicional.

ln(D_cm) = a0 + a1·ln(k_y) + a2·(ln(k_y))² + a3·ln(k_y)·ln(Sa(1.5·T_s))
         + a4·ln(Sa(1.5·T_s)) + a5·M + ε
(coeficientes Bray-Travasarou Tabla 4)
"""
from __future__ import annotations

from math import log, sqrt
from dataclasses import dataclass


@dataclass(frozen=True)
class BrayTravasarouInput:
    k_y: float                       # aceleración crítica (g)
    Sa_1_5Ts_g: float                # Sa a 1.5·T_s (g)
    Mw: float
    Ts_seconds: float = 0.0          # período del talud (0 = rígido)


@dataclass(frozen=True)
class BrayTravasarouResult:
    median_displacement_cm: float
    sigma_lnD: float                 # σ ≈ 0.66
    P_D_gt_5cm: float                # probabilidad de superar 5 cm
    P_D_gt_30cm: float


def evaluate_bray_travasarou(inp: BrayTravasarouInput) -> BrayTravasarouResult:
    """Modelo Bray-Travasarou 2007 (caso Ts > 0.05 s)."""
    if inp.k_y <= 0 or inp.Sa_1_5Ts_g <= 0:
        raise ValueError("k_y y Sa deben ser > 0")

    # Coeficientes Bray-Travasarou 2007 (J. Geotech. ASCE), Eq. (14.5).
    # BUG CORREGIDO (campana de goldens 2026-08-26): la version anterior usaba
    # SIEMPRE el intercepto -0.22 (que es el de la variante RIGIDA, Ts<0.05 s)
    # combinado con el termino flexible +1.50*Ts — para Ts>=0.05 la ecuacion
    # publicada usa -1.10, y la mezcla sobrepredecia el mediano x e^0.88=2.41
    # (caso Austrian Dam: 82.4 cm frente a los 34.2 publicables). La variante
    # rigida ademas no lleva el termino de Ts.
    ln_ky = log(inp.k_y)
    ln_Sa = log(inp.Sa_1_5Ts_g)
    if inp.Ts_seconds < 0.05:
        intercept, ts_term = -0.22, 0.0
    else:
        intercept, ts_term = -1.10, 1.50 * inp.Ts_seconds
    ln_D = (
        intercept
        - 2.83 * ln_ky
        - 0.333 * ln_ky ** 2
        + 0.566 * ln_ky * ln_Sa
        + 3.04 * ln_Sa
        - 0.244 * ln_Sa ** 2
        + ts_term
        + 0.278 * (inp.Mw - 7)
    )
    D_median = float(__import__("math").exp(ln_D))  # cm
    sigma_lnD = 0.66

    # P(D > x) usando lognormal
    import math
    def p_exceed(x_cm: float) -> float:
        if x_cm <= 0:
            return 1.0
        z = (log(x_cm) - ln_D) / sigma_lnD
        # Q(z) = 1 - Φ(z)
        return 0.5 * (1 - math.erf(z / math.sqrt(2)))

    return BrayTravasarouResult(
        median_displacement_cm=D_median,
        sigma_lnD=sigma_lnD,
        P_D_gt_5cm=p_exceed(5.0),
        P_D_gt_30cm=p_exceed(30.0),
    )


__all__ = ["BrayTravasarouInput", "BrayTravasarouResult",
           "evaluate_bray_travasarou"]
