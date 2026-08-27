# -*- coding: utf-8 -*-
"""Golden tests — geotextil en Bishop, anclado a la ficha W00222.

Ficha máster (fuente real):
  W00222 — "Mejora del terreno" (máster UPC), §6 "Drenaje vertical acelerado
  + geotextiles" (pp. 51-59):
    "Geotextiles de refuerzo: aumentan FS en bases de terraplenes sobre
     suelos blandos. Esquema de cálculo: fuerza tensil T_geotextil
     contribuye al momento resistente en análisis Bishop."
  WorkID: W00222. La ficha da la CLÁUSULA (T entra como momento resistente
  en Bishop), no un ejemplo numérico resuelto; la magnitud se ancla con la
  aritmética cerrada de abajo.

¿Por qué NO es tautológico?
  Para φ'=0 el numerador de Bishop es independiente de FS (mα=cosα), así que
  el FS reforzado tiene forma cerrada EXACTA:

      FS_ref - FS_0 = M_geo / (R · M_d) = T·(cy - y_geo) / Σ W·(x_m - cx)

  con brazo vertical (cy - y_geo) porque T es horizontal (reinforcement.py).
  El test calcula el término motor Σ W·(x_m - cx) directamente de la
  discretización (misma vía independiente que test_golden_slip_methods) y el
  momento del geotextil a MANO: T=50 kN/m, brazo = 22-0 = 22 m ->
  M_geo = 1100 kN·m/m. Comparte con la librería únicamente los pesos W
  (mismo alcance declarado en test_golden_slip_methods).

Geometría (idéntica al golden φ'=0 existente): talud 1V:2H, H=10 m,
arcilla c_u=25 kPa, γ=18; círculo cx=15, cy=22, R=24. Geotextil horizontal
en y=0 (base del terraplén — exactamente el caso de la ficha), extensión
x∈[0,30], T_allow=50 kN/m, tau_bond=30 kPa.

Comprobaciones a mano del corte y del anclaje:
  y=0: disc = R²-(cy)² = 576-484 = 92 -> x = 15 ± 9.5917 -> corta el círculo
  en x=5.408 (dentro de [0,30]).  L_anclada = max(5.408, 24.592) = 24.59 m
  -> T_po = 2·30·24.59 = 1476 kN/m >> 50 -> gobierna T_allow = 50 kN/m.
"""
from __future__ import annotations

import pytest

from slopestabilityjm.soil_profile import SoilLayer
from slopestabilityjm.geometry import Slope, CircularSurface
from slopestabilityjm.bishop import BishopAnalysis
from slopestabilityjm.reinforcement import GeosyntheticLayer

_C = 25.0
_CX, _CY, _R = 15.0, 22.0, 24.0
_N = 200
_T_ALLOW = 50.0
_Y_GEO = 0.0


def _slope() -> Slope:
    layer = SoilLayer("Clay", cohesion=_C, friction_angle=0.0,
                      unit_weight=18.0)
    return Slope(profile=[(0, 0), (20, 10), (30, 10)],
                 layers=[layer], layer_depths=[40.0])


def _surface() -> CircularSurface:
    return CircularSurface(cx=_CX, cy=_CY, R=_R)


def _geotextile() -> GeosyntheticLayer:
    return GeosyntheticLayer(y=_Y_GEO, x_left=0.0, x_right=30.0,
                             T_allow=_T_ALLOW, tau_bond=30.0)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestGeotextileBishopW00222:

    def test_reinforcement_moment_hand_value(self):
        """M_geo = T·(cy - y_geo) = 50·22 = 1100 kN·m/m exacto (a mano)."""
        res = BishopAnalysis(_slope(), n_slices=_N,
                             reinforcements=[_geotextile()]
                             ).factor_of_safety(_surface())
        assert res["M_reinforcement_kNm"] == pytest.approx(1100.0, abs=1e-9)

    def test_fs_increase_closed_form(self):
        """FS_ref - FS_0 == T·(cy-y_geo) / Σ W·(x_m-cx)  (cerrada φ'=0).

        El término motor se reconstruye de la discretización (vía
        independiente del bucle de bishop.py; sinα=(x_m-cx)/R exacto).
        """
        slope = _slope()
        fs0 = BishopAnalysis(slope, n_slices=_N).factor_of_safety(
            _surface())["FS"]
        fs_r = BishopAnalysis(slope, n_slices=_N,
                              reinforcements=[_geotextile()]
                              ).factor_of_safety(_surface())["FS"]
        sl = slope.discretize_circle(_surface(), n_slices=_N)
        driving = sum(s["weight_kN"] * (s["x_mid"] - _CX) for s in sl)
        delta_closed = _T_ALLOW * (_CY - _Y_GEO) / driving
        # FS se publica redondeado a 4 decimales -> tolerancia 2.5e-4
        assert fs_r - fs0 == pytest.approx(delta_closed, abs=2.5e-4), (
            f"ΔFS={fs_r - fs0:.6f} != cerrada {delta_closed:.6f} "
            "(W00222: T_geotextil como momento resistente en Bishop)"
        )
        # y el refuerzo efectivamente AUMENTA FS (cláusula de la ficha)
        assert fs_r > fs0

    def test_design_force_governed_by_T_allow(self):
        """Anclaje sobrado: T_po = 2·30·24.59 = 1476 >> 50 -> gobierna
        T_allow (aritmética a mano del docstring)."""
        geo = _geotextile()
        xi = geo.intersect_circle(_surface())
        assert xi == pytest.approx(5.4083, abs=2e-3)  # 15 - sqrt(92)
        assert geo.design_force(xi, sigma_n=50.0) == pytest.approx(50.0)
