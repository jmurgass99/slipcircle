# -*- coding: utf-8 -*-
"""Golden tests — anclaje no tautológico de los métodos de dovelas.

Procedencia de la cita (split honesto WorkID <-> número):
- El **método** y la **cláusula** se anclan a fichas máster:
    * W00254 (Precarga slides, p.39-43): factor de seguridad de talud
      ``FS = Mc/Md``, **método de Bishop simplificado** sobre **círculos de
      deslizamiento**. WorkID: W00254.
    * W01553 (Construction Methods & Project Management): forma cerrada de
      estabilidad de talud ``FS = (c' + σ'·tanφ') / τ``. WorkID: W01553.
    * W01672 (EN 1997-1, EC7): estado límite GEO / **estabilidad global**
      (overall stability, EC7 §11). WorkID: W01672.
- El **valor numérico (magnitud del FS)** NO procede de ninguna ficha: tras
  barrer el `_INDICE_MAESTRO.md` no existe en el corpus un problema de talud
  con la tabla dovela-a-dovela y su FS publicado. Por eso el ancla de magnitud
  es **analítica** (forma cerrada del círculo no drenado), no una cita.

¿Por qué NO es tautológico?
  Para φ'=0 (suelo puramente cohesivo) ambos métodos de dovelas colapsan
  EXACTAMENTE a la forma cerrada del círculo no drenado:

        FS = c · R · L_arco / Σ[ W · (x_m − cx) ]

  La librería calcula el numerador como Σ(c·b/cosα) **dovela a dovela**.
  Este test, en cambio, calcula ``L_arco = R·Δθ`` a partir de la **geometría
  del círculo** (ángulos de entrada/salida), una vía completamente
  independiente del bucle por dovelas.

  Alcance exacto del ancla: en φ'=0 el término motor (ΣW·(xm−cx)) es idéntico
  en ambas vías (sinα=(xm−cx)/R es exacto), por lo que TODO el poder
  discriminante reside en ``L_arco = R·Δθ`` vs ``Σ(b/cosα)``. Es decir, este
  golden guarda la **reducción de la fórmula** (factor mα y geometría de la
  base): atrapa el bug clásico de "olvidar mα" (usar la cuerda 2R·sin(Δθ/2)
  ≈ chord en vez del arco ≈ +6%). NO comprueba de forma independiente los
  pesos W (ambas vías consumen el mismo W). Adicionalmente se exige
  Bishop ≡ Fellenius en φ'=0 (dos reductores escritos de forma independiente
  que coinciden con la reducción analítica conocida).
"""
from __future__ import annotations

import math
from math import asin, cos, degrees

import pytest

from slipcircle.soil_profile import SoilLayer
from slipcircle.geometry import Slope, CircularSurface
from slipcircle.bishop import BishopAnalysis
from slipcircle.fellenius import Slice, fellenius_fs


# --- Caso φ'=0 con geometría/parámetros conocidos -------------------------
# Talud 1V:2H, H=10 m. Círculo elegido para que entre/salga limpiamente
# dentro de [toe_x, crest_x] (evita el clamp espurio del arco).
_C = 25.0          # cohesión [kPa]
_CX, _CY, _R = 15.0, 22.0, 24.0
_N = 200           # dovelas (convergido y estable)


def _phi0_slope() -> Slope:
    layer = SoilLayer("Clay", cohesion=_C, friction_angle=0.0, unit_weight=18.0)
    return Slope(
        profile=[(0, 0), (20, 10), (30, 10)],  # toe -> break -> crest
        layers=[layer], layer_depths=[40.0],
    )


def _surface() -> CircularSurface:
    return CircularSurface(cx=_CX, cy=_CY, R=_R)


def _fellenius_slices(slope: Slope, surface: CircularSurface, n: int):
    """Convierte la discretización de la librería en ``fellenius.Slice``.

    Mapea l = b/cosα (longitud de la base de la dovela)."""
    sl = slope.discretize_circle(surface, n_slices=n)
    return sl, [
        Slice(
            W_kN=s["weight_kN"],
            alpha_deg=degrees(s["alpha_rad"]),
            c_kPa=s["cohesion"],
            phi_deg=degrees(s["phi_rad"]),
            base_length_m=s["width"] / cos(s["alpha_rad"]),
            u_kPa=s["u_kPa"],
        )
        for s in sl
    ]


def _closed_form_fs(slices: list[dict]) -> float:
    """FS analítico del círculo no drenado usando arco GEOMÉTRICO.

    L_arco = R·Δθ donde Δθ se deriva de las posiciones x de entrada/salida
    de las dovelas reales (no del círculo completo), para respetar el clamp
    en toe/crest. Vía independiente del Σ(c·b/cosα) de la librería.
    """
    dx = slices[0]["width"]
    x_lo = slices[0]["x_mid"] - dx / 2.0
    x_hi = slices[-1]["x_mid"] + dx / 2.0
    th_lo = asin((x_lo - _CX) / _R)
    th_hi = asin((x_hi - _CX) / _R)
    L_arc = _R * abs(th_hi - th_lo)
    driving = sum(s["weight_kN"] * (s["x_mid"] - _CX) for s in slices)
    return _C * _R * L_arc / driving


class TestGoldenPhiZeroClosedForm:
    """Ancla de magnitud analítica para φ'=0 (no tautológica)."""

    def test_bishop_matches_undrained_closed_form(self):
        """Bishop (Σ c·b/cosα dovela-a-dovela) == c·R·L_arco/ΣW·(xm−cx)."""
        slope = _phi0_slope()
        fs_bishop = BishopAnalysis(slope, n_slices=_N).factor_of_safety(_surface())["FS"]
        sl = slope.discretize_circle(_surface(), n_slices=_N)
        fs_closed = _closed_form_fs(sl)
        # Tolerancia holgada frente al residuo de discretización (~5e-5) y al
        # redondeo de FS de Bishop a 4 decimales.
        assert fs_bishop == pytest.approx(fs_closed, rel=1e-3), (
            f"Bishop FS={fs_bishop} no coincide con forma cerrada "
            f"FS={fs_closed:.6f} (W00254 Bishop/círculo; ancla analítica)"
        )

    def test_fellenius_matches_undrained_closed_form(self):
        """Fellenius (Σ c·l, l=b/cosα) == forma cerrada analítica."""
        sl, fsl = _fellenius_slices(_phi0_slope(), _surface(), _N)
        fs_fell = fellenius_fs(fsl)
        fs_closed = _closed_form_fs(sl)
        assert fs_fell == pytest.approx(fs_closed, rel=1e-3), (
            f"Fellenius FS={fs_fell:.6f} no coincide con forma cerrada "
            f"FS={fs_closed:.6f}"
        )

    def test_bishop_equals_fellenius_at_phi_zero(self):
        """En φ'=0 ambos reductores (independientes) deben coincidir.

        mα→cosα (Bishop pierde el término en FS) y c·l = c·b/cosα (Fellenius):
        los dos colapsan a la misma reducción analítica conocida.
        """
        _sl, fsl = _fellenius_slices(_phi0_slope(), _surface(), _N)
        fs_fell = fellenius_fs(fsl)
        fs_bishop = BishopAnalysis(_phi0_slope(), n_slices=_N).factor_of_safety(_surface())["FS"]
        assert fs_fell == pytest.approx(fs_bishop, abs=2e-3), (
            f"Fellenius={fs_fell:.6f} vs Bishop={fs_bishop:.6f} divergen en φ'=0"
        )


class TestSliceMethodOrdering:
    """Ordenación física Fellenius <= Bishop para suelo c'-φ' (no φ'=0).

    Cláusula: método de dovelas; el Ordinario (Fellenius) desprecia las
    fuerzas entre dovelas y es CONSERVADOR. WorkID método: W00254.
    """

    def test_fellenius_not_greater_than_bishop(self):
        layer = SoilLayer("CPhi", cohesion=15.0, friction_angle=22.0,
                          unit_weight=19.0)
        slope = Slope(profile=[(0, 0), (20, 10), (30, 10)],
                      layers=[layer], layer_depths=[40.0])
        surf = _surface()
        sl = slope.discretize_circle(surf, n_slices=120)
        fsl = [
            Slice(s["weight_kN"], degrees(s["alpha_rad"]), s["cohesion"],
                  degrees(s["phi_rad"]), s["width"] / cos(s["alpha_rad"]),
                  s["u_kPa"])
            for s in sl
        ]
        fs_fell = fellenius_fs(fsl)
        fs_bishop = BishopAnalysis(slope, n_slices=120).factor_of_safety(surf)["FS"]
        assert fs_fell <= fs_bishop + 1e-9, (
            f"Fellenius={fs_fell:.4f} debe ser <= Bishop={fs_bishop:.4f}"
        )
