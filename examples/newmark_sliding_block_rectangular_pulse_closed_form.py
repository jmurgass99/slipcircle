"""Bloque deslizante de Newmark con pulso rectangular: solución cerrada frente a las dos integraciones del paquete.

Pulso de 0,5 g durante 0,5 s y reposo después, aceleración crítica k_y = 0,2 g. Doble
integración a mano: durante el pulso a_rel = (A − k_y)·g; después el bloque frena a k_y·g
hasta pararse: D = ½·(A − k_y)·g·T²·(A/k_y) = 0,9197 m con g = 9,81 (0,9196 con 9,80665).
Se comparan `evaluate_newmark` (módulo `newmark`) y `PseudoStaticAnalysis.newmark_integration`,
y se enseña que el desplazamiento no depende de cuánto reposo se añada al final (el defecto
que tenía el módulo hasta la campaña de goldens de agosto de 2026: el bloque «planeaba»).

Ejecutar:  python examples/newmark_sliding_block_rectangular_pulse_closed_form.py
"""
import numpy as np

from slipcircle import PseudoStaticAnalysis
from slipcircle.newmark import NewmarkInput, evaluate_newmark

A, KY, T, DT = 0.5, 0.2, 0.5, 0.001


def main() -> None:
    for g in (9.81, 9.80665):
        print(f"Forma cerrada con g = {g}: D = {0.5 * (A - KY) * g * T ** 2 * (A / KY):.5f} m")
    for tail_s in (1.0, 2.0, 4.0):
        acc = np.concatenate([np.full(int(T / DT), A), np.zeros(int(tail_s / DT))])
        r = evaluate_newmark(NewmarkInput(accelerogram_g=acc.tolist(), dt_s=DT, a_y_g=KY))
        t = np.arange(acc.size) * DT
        r2 = PseudoStaticAnalysis.newmark_integration(t, acc, ky_g=KY)
        peak = r2.get("peak_displacement_m", r2.get("displacement_m", max(r2["displacement"]) if "displacement" in r2 else float("nan")))
        print(f"cola de reposo {tail_s:.0f} s: evaluate_newmark D = {r.permanent_displacement_m:.5f} m "
              f"({r.n_pulses_exceeding} pulso, I_a = {r.arias_intensity_m_s:.3f} m/s, {r.severity}); "
              f"newmark_integration D = {peak:.5f} m")


if __name__ == "__main__":
    main()
