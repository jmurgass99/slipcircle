"""
Newmark sliding block analysis — desplazamiento permanente acumulado en
un talud bajo acelerograma sísmico, integrando los pulsos que exceden la
aceleración crítica a_y.

Referencias:
    Newmark, N.M. (1965) — Effects of earthquakes on dams and embankments.
    Jibson, R.W. (2007) — Regression models for estimating coseismic
    landslide displacement.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class NewmarkInput:
    accelerogram_g: Sequence[float]    # registro del suelo (g)
    dt_s: float                         # paso de tiempo
    a_y_g: float                        # aceleración crítica del talud (g)


@dataclass(frozen=True)
class NewmarkResult:
    permanent_displacement_m: float
    n_pulses_exceeding: int
    arias_intensity_m_s: float
    severity: str  # 'low'<1cm, 'medium'<10cm, 'high'>=10cm


def evaluate_newmark(inp: NewmarkInput) -> NewmarkResult:
    if inp.dt_s <= 0 or inp.a_y_g <= 0:
        raise ValueError("dt_s y a_y_g deben ser > 0")
    ag = np.asarray(inp.accelerogram_g, dtype=float)
    if ag.size == 0:
        raise ValueError("accelerogram_g vacío")

    g = 9.80665
    a_y = inp.a_y_g * g  # m/s²
    accel = ag * g

    # Integrar pulsos
    velocity = 0.0
    displacement = 0.0
    n_exc = 0
    in_pulse = False
    for a in accel:
        if abs(a) > a_y:
            # Aceleración relativa hacia abajo (signo convención: a > a_y desliza)
            a_rel = (abs(a) - a_y) * np.sign(a)
            velocity += a_rel * inp.dt_s
            if not in_pulse:
                in_pulse = True
                n_exc += 1
        if velocity > 0:
            displacement += velocity * inp.dt_s
        else:
            velocity = 0.0  # no se mueve hacia atrás
            in_pulse = False

    # Arias intensity (m/s)
    Ia = (np.pi / (2 * g)) * np.sum(accel ** 2) * inp.dt_s

    if displacement < 0.01:
        sev = "low"
    elif displacement < 0.10:
        sev = "medium"
    else:
        sev = "high"

    return NewmarkResult(
        permanent_displacement_m=float(displacement),
        n_pulses_exceeding=n_exc,
        arias_intensity_m_s=float(Ia),
        severity=sev,
    )


__all__ = ["NewmarkInput", "NewmarkResult", "evaluate_newmark"]
