"""Desplazamiento permanente de Newmark para un acelerograma sintético y varias aceleraciones críticas.

Registro sintético de 20 s (suma de armónicos de 1 a 4 Hz con envolvente de Saragoni-Hart,
PGA escalada a 0,35 g): `evaluate_newmark` integra el bloque deslizante para k_y de 0,05 a
0,30 g y clasifica la severidad (< 1 cm baja, < 10 cm media, ≥ 10 cm alta). La aceleración
crítica del talud se obtiene del análisis pseudo-estático como el k_h que lleva FS a 1,0.

Ejecutar:  python examples/newmark_displacement_from_synthetic_accelerogram_ky_sweep.py
"""
import numpy as np

from slipcircle.newmark import NewmarkInput, evaluate_newmark

DT, DUR, PGA = 0.01, 20.0, 0.35


def synthetic_record(seed: int = 4) -> np.ndarray:
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, DUR, DT)
    env = (t / 4.0) ** 2 * np.exp(-t / 4.0)
    sig = sum(a * np.sin(2 * np.pi * f * t + p) for a, f, p in zip(rng.uniform(0.3, 1.0, 8), rng.uniform(1.0, 4.0, 8), rng.uniform(0, 2 * np.pi, 8)))
    acc = env * sig
    return PGA * acc / np.max(np.abs(acc))


def main() -> None:
    acc = synthetic_record()
    print(f"Registro: {acc.size} puntos, PGA = {np.max(np.abs(acc)):.3f} g, Δt = {DT} s\n")
    print(f"{'k_y [g]':>8} {'D [cm]':>8} {'pulsos':>7} {'I_a [m/s]':>10}  severidad")
    for ky in (0.05, 0.10, 0.15, 0.20, 0.25, 0.30):
        r = evaluate_newmark(NewmarkInput(accelerogram_g=acc.tolist(), dt_s=DT, a_y_g=ky))
        print(f"{ky:8.2f} {100 * r.permanent_displacement_m:8.2f} {r.n_pulses_exceeding:7d} {r.arias_intensity_m_s:10.3f}  {r.severity}")
    print("\nSolo cuentan los pulsos de signo ladera abajo que superan k_y; por encima de la PGA no hay deslizamiento.")


if __name__ == "__main__":
    main()
