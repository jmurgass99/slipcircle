# -*- coding: utf-8 -*-
"""Resolución de (FS, p) en los métodos con fuerzas entre dovelas paralelas (Spencer)
o de inclinación prescrita (Morgenstern-Price, p = λ).

Esquema clásico de Spencer (1967): para cada valor del parámetro p se obtiene
FS_f(p) con el equilibrio de fuerzas y FS_m(p) con el de momentos, cada uno como
raíz unidimensional en FS; la solución es el p en que FS_f = FS_m. Es el
procedimiento de los textos (Duncan & Wright, cap. 6) y es mucho más robusto que
alternar raíces en FS y en p, que es lo que hacía este paquete hasta 2026-10-31:
el resultado dependía del número de dovelas (θ se quedaba en 0 sin satisfacer el
equilibrio de momentos, o se disparaba a FS = 8,7 con θ = −46°).

La función Q_i(FS, p) de cada dovela tiene un polo en
FS = −tan(α_i − θ_i)·tan φ_i (donde 1 + tan(α−θ)·tan φ/FS se anula). Por eso el
extremo inferior del bracket de FS se sitúa, para cada p, por encima del mayor de
esos polos (``fs_lower(p)``): un brentq que lo incluyera devolvería una raíz falsa.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import brentq

Residuals = Callable[[float, float], tuple[float, float]]


@dataclass(frozen=True)
class ParallelSolution:
    FS: float
    p: float
    FS_force: float | None
    FS_moment: float | None
    converged: bool
    n_scanned: int
    n_valid: int


def solve_parallel(residuals: Residuals, p_grid: Sequence[float],
                   fs_lower: Callable[[float], float], fs_hi: float = 10.0,
                   tol: float = 1e-6) -> ParallelSolution:
    """Busca p tal que la raíz en FS del residuo de fuerzas coincide con la del de momentos.

    ``residuals(FS, p)`` devuelve ``(sumF, sumM)``. ``fs_lower(p)`` es el mayor
    polo de Q en FS para ese p (0 si no hay). Si en ningún tramo de la malla hay
    cambio de signo de FS_f − FS_m, se devuelve el p de menor |FS_f − FS_m| con
    ``converged=False``: el llamador debe mirar la bandera, no solo el FS.
    """
    def root(idx: int, p: float) -> float | None:
        lo = max(0.3, 1.1 * float(fs_lower(p)))
        if lo >= fs_hi:
            return None

        def f(fs: float) -> float:
            return residuals(fs, p)[idx]

        try:
            flo, fhi = f(lo), f(fs_hi)
        except (ZeroDivisionError, OverflowError, ValueError):
            return None
        if not (np.isfinite(flo) and np.isfinite(fhi)) or flo * fhi > 0:
            return None
        return float(brentq(f, lo, fs_hi, xtol=tol))

    rows: list[tuple[float, float | None, float | None]] = []
    for p in p_grid:
        rows.append((float(p), root(0, float(p)), root(1, float(p))))
    valid = [(p, ff, fm) for p, ff, fm in rows if ff is not None and fm is not None]
    if not valid:
        return ParallelSolution(float("nan"), float("nan"), None, None, False, len(rows), 0)

    def diff(p: float) -> float:
        ff, fm = root(0, p), root(1, p)
        if ff is None or fm is None:
            raise ValueError("sin raíz en FS para este p")
        return ff - fm

    for (p0, f0, m0), (p1, f1, m1) in zip(valid, valid[1:]):
        d0, d1 = f0 - m0, f1 - m1
        if d0 == 0.0:
            return ParallelSolution(f0, p0, f0, m0, True, len(rows), len(valid))
        if d0 * d1 < 0.0:
            try:
                p_star = float(brentq(diff, p0, p1, xtol=1e-7))
            except ValueError:
                p_star = p0 + (p1 - p0) * d0 / (d0 - d1)
            ff, fm = root(0, p_star), root(1, p_star)
            if ff is not None and fm is not None:
                conv = abs(ff - fm) < 1e-4 * max(ff, 1.0)
                return ParallelSolution(0.5 * (ff + fm), p_star, ff, fm, conv, len(rows), len(valid))
    p, ff, fm = min(valid, key=lambda r: abs(r[1] - r[2]))
    return ParallelSolution(0.5 * (ff + fm), p, ff, fm, False, len(rows), len(valid))


__all__ = ["ParallelSolution", "solve_parallel"]
