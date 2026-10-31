# -*- coding: utf-8 -*-
"""Golden tests — Bray & Travasarou (2007) contra la ecuación PUBLICADA.

Fuente literal verificada (2026-08-26):
  Bray, J.D. (2007) "Simplified seismic slope displacement procedures",
  cap. 14 de *Earthquake Geotechnical Engineering* (Springer), §5.5,
  Eq. (14.5), pp. 344 — idéntica a Bray & Travasarou (2007) JGGE 133(4):

    ln(D) = -1.10 - 2.83·ln(k_y) - 0.333·(ln k_y)²
            + 0.566·ln(k_y)·ln(Sa(1.5Ts)) + 3.04·ln(Sa(1.5Ts))
            - 0.244·(ln Sa(1.5Ts))² + 1.5·Ts + 0.278·(M-7) ± ε,  σ = 0.66

  y cita textual: "the first term of Eq. (14.5) should be replaced with
  −0.22 when Ts < 0.05 s" (caso bloque casi rígido, Sa -> PGA).

Fichas del corpus: NINGUNA menciona Bray-Travasarou (barrido de
`_EXTRACCION_DOCUMENTOS/fichas/` 2026-08-26). Ancla = ecuación publicada.

BUG DE PRODUCCIÓN documentado (NO corregido, mandato de campaña):
  bray_travasarou.py:41 usa el intercepto RÍGIDO (-0.22) y a la vez
  bray_travasarou.py:47 suma el término FLEXIBLE (+1.50·Ts), con docstring
  "(caso Ts > 0.05 s)". Para Ts >= 0.05 s el mediano queda multiplicado por
  exp(-0.22-(-1.10)) = exp(0.88) = 2.411 respecto al valor publicado.
  Para Ts = 0 (bloque rígido) la implementación SÍ es correcta.
"""
from __future__ import annotations

import math

import pytest

from slipcircle.bray_travasarou import (
    BrayTravasarouInput,
    evaluate_bray_travasarou,
)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestRigidBlockPublishedEquation:
    """Caso rígido Ts=0 — la implementación coincide con la variante
    publicada (-0.22, sin término en Ts porque Ts=0)."""

    def test_rigid_case_hand_value(self):
        """k_y=0.10, PGA=Sa=0.40 g, M=7.0, Ts=0.

        Aritmética a mano (Eq. 14.5 con -0.22, coeficientes del paper):
          ln k_y = -2.302585 ; ln Sa = -0.916291
          -0.22 + 6.516316 - 1.765532 + 1.194173 - 2.785525 - 0.204860
          = 2.734572  ->  D = e^2.734572 = 15.403 cm
        """
        res = evaluate_bray_travasarou(
            BrayTravasarouInput(k_y=0.10, Sa_1_5Ts_g=0.40, Mw=7.0,
                                Ts_seconds=0.0))
        assert res.median_displacement_cm == pytest.approx(15.403, rel=1e-3)
        assert res.sigma_lnD == pytest.approx(0.66)

    def test_probability_lognormal_hand_value(self):
        """P(D > 5 cm) con mediana 15.403 cm y σ=0.66.

        z = (ln 5 - 2.734572)/0.66 = -1.7047; tabla normal estándar:
        P = Φ(1.7047) = 0.9559 (vía independiente: tabla, no erf del módulo).
        """
        res = evaluate_bray_travasarou(
            BrayTravasarouInput(k_y=0.10, Sa_1_5Ts_g=0.40, Mw=7.0,
                                Ts_seconds=0.0))
        assert res.P_D_gt_5cm == pytest.approx(0.9559, abs=2e-3)

    def test_half_to_double_range_property(self):
        """Propiedad publicada: exp(σ)=exp(0.66)≈1.93≈2 -> el rango 16-84%
        es ~ mitad/doble de la mediana (Bray 2007, p. 344)."""
        assert math.exp(0.66) == pytest.approx(1.935, abs=5e-3)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestFlexibleCaseInterceptBug:
    """Caso flexible (Ts >= 0.05 s) — el intercepto implementado NO es el
    publicado. Datos del caso de validación 'Austrian Dam' (Loma Prieta),
    Tabla 14.1 del propio capítulo: k_y=0.14, Ts=0.33 s, Sa(1.5Ts)=0.94 g;
    rango B&T publicado 20-70 cm (Dmax observado 50 cm).

    Aritmética a mano con Eq. (14.5), M=6.93 (Loma Prieta):
      ln(D) = -1.10 + 5.5641 - 1.2872 + 0.0689 - 0.1881 - 0.0009
              + 0.4950 - 0.0195 = 3.5322  ->  D_mediana = 34.2 cm  (en rango)
    La implementación devuelve 82.45 cm = 2.411 × 34.2 (fuera del rango
    publicado 20-70 cm): factor exp(0.88) del intercepto erróneo.
    """

    def test_austrian_dam_median_matches_published_equation(self):
        """BUG CORREGIDO (antes xfail estricto): la version previa usaba el
        intercepto -0.22 (rigido) con el termino flexible +1.50*Ts y daba
        82.4 cm; la Eq. (14.5) publicada (intercepto -1.10 para Ts>=0.05 s)
        da 34.2 cm para el caso Austrian Dam — dentro del rango 20-70 cm de
        la Tabla 14.1 del propio paper."""
        res = evaluate_bray_travasarou(
            BrayTravasarouInput(k_y=0.14, Sa_1_5Ts_g=0.94, Mw=6.93,
                                Ts_seconds=0.33))
        assert res.median_displacement_cm == pytest.approx(34.2, rel=0.02)

    def test_flexible_case_matches_published_equation_exactly(self):
        """Guardia del fix (antes 'evidencia del bug', que afirmaba el factor
        exp(0.88)): el mediano implementado coincide 1:1 con la Eq. (14.5)
        transcrita a mano para el caso flexible."""
        res = evaluate_bray_travasarou(
            BrayTravasarouInput(k_y=0.14, Sa_1_5Ts_g=0.94, Mw=6.93,
                                Ts_seconds=0.33))
        ln_ky = math.log(0.14)
        ln_sa = math.log(0.94)
        ln_D_published = (
            -1.10 - 2.83 * ln_ky - 0.333 * ln_ky ** 2
            + 0.566 * ln_ky * ln_sa + 3.04 * ln_sa - 0.244 * ln_sa ** 2
            + 1.50 * 0.33 + 0.278 * (6.93 - 7.0)
        )
        ratio = res.median_displacement_cm / math.exp(ln_D_published)
        assert ratio == pytest.approx(1.0, rel=1e-9)
