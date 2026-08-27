# -*- coding: utf-8 -*-
"""Golden tests — bloque deslizante de Newmark con pulso rectangular.

Procedencia (split honesto WorkID <-> número):
- **Método**: Newmark, N.M. (1965) "Effects of earthquakes on dams and
  embankments" (referencia citada en ``slopestabilityjm/newmark.py`` y
  descrita como "Newmark rigid sliding block procedure" en Bray, J.D. (2007)
  *Simplified seismic slope displacement procedures*, cap. 14, §4.4).
- **Fichas del corpus**: NO existe ficha con ejemplo resuelto del bloque
  deslizante. Todas las apariciones de "Newmark" en las fichas del máster
  (W00018, W00019, W01539, W00882, ...) son la **integración temporal
  Newmark-β** de dinámica estructural — un concepto distinto que NO ancla
  este módulo. Se deja constancia para no confundir anclas.
- **Magnitud**: solución CERRADA del pulso rectangular, derivada a mano
  (doble integración), vía completamente independiente del bucle de
  integración de la librería.

Solución cerrada (pulso rectangular a = A·g durante T, luego reposo):
  Fase 1 (0..T):  a_rel = (A-k_y)·g constante
                  v(T) = (A-k_y)·g·T ;  d(T) = ½·(A-k_y)·g·T²
  Fase 2 (t>T):   el bloque desacelera a ritmo k_y·g (fricción) hasta v=0
                  d_extra = v(T)²/(2·k_y·g)
  TOTAL:          D = ½·(A-k_y)·g·T²·(A/k_y)

Con A=0.5 g, k_y=0.2 g, T=0.5 s:
  v(T) = 0.3·g·0.5 = 1.4715 m/s   (g=9.81)
  D_pulso = ½·0.3·9.81·0.25 = 0.36788 m  (g=9.81) / 0.36775 m (g=9.80665)
  D_total = D_pulso·(0.5/0.2) = 0.9196875 m (g=9.81)
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from slopestabilityjm.newmark import NewmarkInput, evaluate_newmark
from slopestabilityjm.seismic import PseudoStaticAnalysis

# Pulso rectangular de referencia
_A = 0.5      # amplitud del pulso [g]
_KY = 0.2     # aceleración crítica [g]
_T = 0.5      # duración del pulso [s]
_DT = 0.001   # paso [s]
_N_PULSE = int(_T / _DT)   # 500 muestras
_G_NEWMARK = 9.80665       # g usado por newmark.py
_G_SEISMIC = 9.81          # g usado por seismic.py (inconsistencia menor)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestNewmarkIntegrationClosedForm:
    """seismic.py::newmark_integration — implementación COMPLETA del bloque
    (desacelera a k_y·g bajo el umbral). Ancla: solución cerrada total."""

    def test_rectangular_pulse_total_displacement(self):
        """D_total == ½(A-k_y)gT²·(A/k_y): pulso + frenado por fricción.

        Cerrada a mano: ½·0.3·9.81·0.25·2.5 = 0.9196875 m. El registro
        lleva 2.0 s de cola en reposo (> t_stop = v_max/(k_y g) = 0.75 s)
        para que el bloque frene por completo dentro del registro.
        """
        t = np.arange(0.0, 2.5 + _DT, _DT)
        accel_g = np.where(t <= _T, _A, 0.0)
        res = PseudoStaticAnalysis.newmark_integration(t, accel_g, _KY)
        D_closed = 0.5 * (_A - _KY) * _G_SEISMIC * _T ** 2 * (_A / _KY)
        assert res["peak_displacement_m"] == pytest.approx(D_closed, rel=5e-3), (
            f"D={res['peak_displacement_m']:.6f} m != cerrada {D_closed:.6f} m "
            "(Newmark 1965, pulso rectangular)"
        )

    def test_velocity_returns_to_zero_after_pulse(self):
        """Tras el frenado el bloque queda parado (v=0): desplazamiento
        permanente, no creciente con la longitud del registro."""
        t = np.arange(0.0, 2.5 + _DT, _DT)
        accel_g = np.where(t <= _T, _A, 0.0)
        res = PseudoStaticAnalysis.newmark_integration(t, accel_g, _KY)
        assert res["velocity_m_s"][-1] == pytest.approx(0.0, abs=1e-9)
        # el desplazamiento final coincide con el pico (ya no crece)
        assert res["displacement_m"][-1] == pytest.approx(
            res["peak_displacement_m"], rel=1e-12)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestEvaluateNewmarkPulsePartials:
    """newmark.py::evaluate_newmark — partes CORRECTAS de la implementación,
    ancladas a la solución cerrada del tramo de pulso y a la intensidad de
    Arias cerrada."""

    def test_pulse_only_displacement(self):
        """Registro que TERMINA al acabar el pulso: la doble integral del
        exceso es ½(A-k_y)gT² (residuo Euler (1+1/N) con N=500 -> +0.2%).

        Cerrada a mano (g=9.80665): ½·0.3·9.80665·0.25 = 0.367749 m.
        """
        res = evaluate_newmark(NewmarkInput([_A] * _N_PULSE, _DT, _KY))
        d_closed = 0.5 * (_A - _KY) * _G_NEWMARK * _T ** 2
        assert res.permanent_displacement_m == pytest.approx(d_closed, rel=3e-3)
        assert res.n_pulses_exceeding == 1
        assert res.severity == "high"  # 0.368 m >= 0.10 m

    def test_arias_intensity_closed_form(self):
        """I_a = π·a²·T/(2g) para pulso rectangular.

        Con a = 0.5·g, T = 0.5 s: I_a = π·0.25·g·T/2 = π·0.25·9.80665·0.25
        = 1.92553 m/s (la suma discreta es exacta para señal constante).
        """
        res = evaluate_newmark(NewmarkInput([_A] * _N_PULSE, _DT, _KY))
        Ia_closed = math.pi * (_A ** 2) * _G_NEWMARK * _T / 2.0
        assert res.arias_intensity_m_s == pytest.approx(Ia_closed, rel=1e-6)

    def test_below_yield_no_displacement(self):
        """a < a_y en todo el registro => D = 0 (cerrada trivial)."""
        res = evaluate_newmark(NewmarkInput([0.15] * 1000, _DT, _KY))
        assert res.permanent_displacement_m == 0.0
        assert res.n_pulses_exceeding == 0


@pytest.mark.unit
@pytest.mark.safety_critical
class TestEvaluateNewmarkMissingDeceleration:
    """BUG DE PRODUCCIÓN documentado (NO corregido, mandato de campaña).

    newmark.py:50-61 — cuando |a| <= a_y y v > 0 la implementación NO aplica
    la desaceleración por fricción k_y·g: el bloque "planea" a velocidad
    constante y el desplazamiento crece LINEALMENTE con la longitud del
    registro (no acotado).

    Aritmética (pulso A=0.5g, k_y=0.2g, T=0.5 s, dt=1 ms):
      física (cerrada): D_total = ½·0.3·9.80665·0.25·(0.5/0.2) = 0.9194 m
      implementación:   D(cola 2 s) = 0.3685 + 1.4710·2.0 = 3.3105 m  (3.6x)
                        D(cola 4 s) = 0.3685 + 1.4710·4.0 = 6.2525 m  (6.8x)
    La implementación correcta del MISMO modelo está en
    seismic.py::newmark_integration (ver TestNewmarkIntegrationClosedForm).
    """

    def test_total_displacement_with_quiet_tail(self):
        """BUG CORREGIDO (antes xfail estricto): la version previa de
        newmark.py no aplicaba la desaceleracion k_y*g con |a|<=a_y y v>0.
        Ahora la cerrada total D = 1/2(A-k_y)gT^2*(A/k_y) se reproduce con
        cola en reposo de 2 s (fisica de Newmark 1965)."""
        record = [_A] * _N_PULSE + [0.0] * 2000
        res = evaluate_newmark(NewmarkInput(record, _DT, _KY))
        D_closed = 0.5 * (_A - _KY) * _G_NEWMARK * _T ** 2 * (_A / _KY)
        assert res.permanent_displacement_m == pytest.approx(D_closed, rel=0.02)

    def test_displacement_is_independent_of_quiet_tail_length(self):
        """Guardia del fix (antes 'evidencia del bug', que afirmaba el
        crecimiento lineal con la cola): la fisica exige que, agotada la
        velocidad, alargar la cola en reposo NO cambie el desplazamiento."""
        d2 = evaluate_newmark(
            NewmarkInput([_A] * _N_PULSE + [0.0] * 2000, _DT, _KY)
        ).permanent_displacement_m
        d4 = evaluate_newmark(
            NewmarkInput([_A] * _N_PULSE + [0.0] * 4000, _DT, _KY)
        ).permanent_displacement_m
        assert d4 == pytest.approx(d2, rel=1e-12)
