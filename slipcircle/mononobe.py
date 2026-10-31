"""
Análisis pseudo-estático Mononobe-Okabe para taludes y muros sísmicos.

FS_seismic = (c·L + N'·tan(φ)) / (W·sin(α+ψ) / cos(ψ))
donde ψ = atan(kh / (1 - kv))
"""
from __future__ import annotations

from dataclasses import dataclass
from math import atan, cos, degrees, radians, sin, tan


@dataclass(frozen=True)
class MononobeInput:
    slope_height_m: float
    slope_angle_deg: float
    gamma_kN_m3: float
    c_kPa: float
    phi_deg: float
    kh: float = 0.15
    kv: float = 0.0


@dataclass(frozen=True)
class MononobeResult:
    FS_static: float
    FS_seismic: float
    psi_deg: float
    is_safe: bool


def evaluate_mononobe(inp: MononobeInput, FS_min: float = 1.1
                        ) -> MononobeResult:
    H = inp.slope_height_m
    alpha = radians(inp.slope_angle_deg)
    phi = radians(inp.phi_deg)
    gamma = inp.gamma_kN_m3
    c = inp.c_kPa

    # Cuña infinita (Bishop simplificado en 2D)
    # FS = (c·L + (γ·H·cos²α - u·L)·tan φ) / (γ·H·sin α·cos α)
    L = 1.0 / cos(alpha)  # longitud unitaria de la base
    W = gamma * H * 1.0   # peso por metro de talud unitario
    N = W * cos(alpha)
    T_drive = W * sin(alpha)
    T_resist = c * L + N * tan(phi)
    FS_static = T_resist / max(T_drive, 1e-9)

    # Pseudo-estático Mononobe
    psi = atan(inp.kh / max(1.0 - inp.kv, 1e-9))
    # Fuerza horizontal kh·W rota el polígono de fuerzas
    Wnew = W * (1.0 - inp.kv) / cos(psi)
    T_drive_seismic = Wnew * sin(alpha + psi)
    N_seismic = Wnew * cos(alpha + psi)
    T_resist_seismic = c * L + N_seismic * tan(phi)
    FS_seismic = T_resist_seismic / max(T_drive_seismic, 1e-9)

    return MononobeResult(
        FS_static=FS_static, FS_seismic=FS_seismic,
        psi_deg=float(degrees(psi)),
        is_safe=FS_seismic >= FS_min,
    )


__all__ = ["MononobeInput", "MononobeResult", "evaluate_mononobe"]
