"""Talud infinito seco con superficie de rotura plana: FS = tan φ / tan β exacto (y con cohesión).

Ladera recta de pendiente β = 15° con plano de rotura paralelo a z = 3 m de profundidad
vertical, φ' = 30°, γ = 20 kN/m³. Es el único caso que `BishopAnalysis` no puede representar
(solo acepta círculos): se discretiza la polilínea con `slices_from_polyline` y se reduce con
el método Ordinario (`fellenius_fs`). Todas las dovelas son idénticas, así que el resultado es
exacto e independiente del número de dovelas: 1 + 2/√3 = 2,1547 sin cohesión y
(5 + 2√3)/3 = 2,8214 con c' = 10 kPa (anclaje de `tests/test_golden_infinite_slope.py`).

Ejecutar:  python examples/infinite_slope_planar_surface_fellenius_polyline.py
"""
from math import cos, degrees, radians, sqrt, tan

from slipcircle import Slope, SoilLayer, slices_from_polyline
from slipcircle.fellenius import Slice, fellenius_fs

BETA, PHI, GAMMA, Z = 15.0, 30.0, 20.0, 3.0


def fs_infinite(c_kPa: float, n_slices: int) -> float:
    tb = tan(radians(BETA))
    slope = Slope(profile=[(0.0, 0.0), (100.0, 100.0 * tb)],
                  layers=[SoilLayer("ladera", cohesion=c_kPa, friction_angle=PHI, unit_weight=GAMMA)],
                  layer_depths=[50.0])
    polyline = [(20.0, 20.0 * tb - Z), (80.0, 80.0 * tb - Z)]
    sl = slices_from_polyline(slope, polyline, n_slices=n_slices)
    return fellenius_fs([Slice(s["weight_kN"], degrees(s["alpha_rad"]), s["cohesion"], degrees(s["phi_rad"]),
                               s["width"] / cos(s["alpha_rad"]), s["u_kPa"]) for s in sl])


def main() -> None:
    print(f"β = {BETA}°, φ' = {PHI}°, γ = {GAMMA} kN/m³, z = {Z} m\n")
    print(f"{'c [kPa]':>8} {'dovelas':>8} {'FS':>10} {'forma cerrada':>14}")
    for c, ref in ((0.0, 1 + 2 / sqrt(3)), (10.0, (5 + 2 * sqrt(3)) / 3)):
        for n in (5, 40, 400):
            print(f"{c:8.0f} {n:8d} {fs_infinite(c, n):10.6f} {ref:14.6f}")
    print("\nSin cohesión FS = tan φ / tan β: no depende de γ ni de z; con cohesión sí (término c/(γ·z·cos β·sin β)).")


if __name__ == "__main__":
    main()
