"""Morgenstern-Price: efecto de la función f(x) entre dovelas (seno, constante, trapezoidal) y Spencer.

Mismo talud y círculo crítico de Duncan y Wright (c' = 10 kPa, φ' = 20°, γ = 18, H = 10 m,
β = 30°). Morgenstern-Price resuelve FS y λ para una f(x) dada; con f constante coincide con
Spencer (ángulo θ constante). Se comparan las tres funciones del paquete: el FS cambia en la
cuarta cifra decimal y lo que cambia es λ, como documenta la literatura para superficies
circulares. Hasta la 1.3.1 el medio seno daba 1,3674 porque no se imponía el equilibrio
vertical de las dovelas.

Ejecutar:  python examples/morgenstern_price_interslice_function_comparison.py
"""
from math import radians, tan

from slipcircle import (
    CircularSurface,
    MorgensternPriceAnalysis,
    Slope,
    SoilLayer,
    SpencerAnalysis,
    constant_f,
    half_sine,
    trapezoidal,
)

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    surf = CircularSurface(3.189, 16.909, 18.091)
    print(f"{'método / f(x)':<28} {'FS':>8} {'λ o θ':>8}  convergió")
    for name, f in (("M-P, medio seno", half_sine), ("M-P, constante", constant_f), ("M-P, trapezoidal", trapezoidal)):
        r = MorgensternPriceAnalysis(slope, n_slices=40, f=f).factor_of_safety(surf)
        print(f"{name:<28} {r['FS']:8.4f} {r['lambda']:8.4f}  {r['converged']}")
    r = SpencerAnalysis(slope, n_slices=40).factor_of_safety(surf)
    theta = r.get("theta_deg", r.get("theta", float("nan")))
    print(f"{'Spencer':<28} {r['FS']:8.4f} {theta:8.3f}  {r['converged']}")


if __name__ == "__main__":
    main()
