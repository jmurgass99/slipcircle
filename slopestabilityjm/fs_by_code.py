"""
Factores de seguridad mínimos para análisis de estabilidad de taludes
por código.
"""
from __future__ import annotations


# Mínimos típicos según norma
FS_MIN_TABLE = {
    # (code, scenario) -> FS_min
    ("EU", "static_persistent"):  1.0,    # γ_R,e via EN 1997
    ("EU", "seismic"):             1.0,    # E_d/R_d ≤ 1
    ("SG", "static_persistent"):  1.0,
    ("US", "static_persistent"):  1.5,    # USACE / FHWA típico
    ("US", "seismic"):             1.1,
    ("CO", "static_persistent"):  1.5,    # NSR-10 §H.5
    ("CO", "seismic"):             1.05,
    ("AU", "static_persistent"):  1.5,    # AS 4678 / GNGS
    ("AU", "seismic"):             1.1,
}


def fs_minimum(code: str, scenario: str = "static_persistent") -> float:
    """Devuelve FS_min según el código (1.0 = EC7 LRFD style)."""
    key = (code.upper(), scenario)
    if key not in FS_MIN_TABLE:
        raise ValueError(f"sin FS_min para {key}")
    return FS_MIN_TABLE[key]


def is_safe(FS: float, code: str, scenario: str = "static_persistent"
             ) -> bool:
    return FS >= fs_minimum(code, scenario)


__all__ = ["FS_MIN_TABLE", "fs_minimum", "is_safe"]
