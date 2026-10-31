# -*- coding: utf-8 -*-
"""Metodo Ordinario de Dovelas (Fellenius / Sueco).

El metodo de equilibrio limite mas simple: desprecia las fuerzas entre dovelas,
por lo que es CONSERVADOR (FS_Fellenius <= FS_Bishop tipicamente). Util como cota
inferior y base de comparacion. Completa la progresion Ordinario -> Bishop ->
Janbu -> Spencer -> Morgenstern-Price.

  FS = Σ[ c'·l + (W·cosα − u·l)·tanφ' ] / Σ[ W·sinα ]

(superficie circular; α = inclinacion de la base de la dovela; u = presion de
poros; l = longitud de la base de la dovela.)
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class Slice:
    W_kN: float           # peso de la dovela (por metro de ancho)
    alpha_deg: float      # inclinacion de la base (+ hacia abajo del talud)
    c_kPa: float          # cohesion efectiva
    phi_deg: float        # angulo de friccion efectivo
    base_length_m: float  # longitud de la base de la dovela
    u_kPa: float = 0.0    # presion de poros en la base


def fellenius_fs(slices: List[Slice], seismic_kh: float = 0.0) -> float:
    """Factor de seguridad por el metodo Ordinario de Dovelas.

    seismic_kh: coeficiente sismico pseudo-estatico horizontal (anade W·kh·... a la
    fuerza motriz, opcional). Devuelve FS.
    """
    if not slices:
        raise ValueError("se requiere al menos una dovela")
    resist = 0.0
    driving = 0.0
    for s in slices:
        a = math.radians(s.alpha_deg)
        phi = math.radians(s.phi_deg)
        N_eff = s.W_kN * math.cos(a) - s.u_kPa * s.base_length_m
        resist += s.c_kPa * s.base_length_m + max(N_eff, 0.0) * math.tan(phi)
        driving += s.W_kN * math.sin(a)
        # Componente sismica pseudo-estatica (fuerza horizontal kh·W actuando
        # con brazo ~ cos? aproximacion estandar: anade kh·W·cosα a la motriz).
        driving += seismic_kh * s.W_kN * math.cos(a)
    if driving <= 0:
        return float("inf")
    return resist / driving


__all__ = ["Slice", "fellenius_fs"]
