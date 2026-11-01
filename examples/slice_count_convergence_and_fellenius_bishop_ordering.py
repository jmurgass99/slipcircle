"""Convergencia con el número de dovelas y la cota conservadora del método Ordinario (Fellenius ≤ Bishop).

Talud de Duncan y Wright sobre su círculo crítico: FS de Bishop con 10 a 320 dovelas, y
Fellenius con 160 dovelas construidas con la misma discretización. En suelo c'-φ' el
Ordinario desprecia las fuerzas entre dovelas y queda del lado de la seguridad; en φ' = 0
coincidirían.

Ejecutar:  python examples/slice_count_convergence_and_fellenius_bishop_ordering.py
"""
from math import cos, degrees, radians, tan

from slipcircle import BishopAnalysis, CircularSurface, Slope, SoilLayer
from slipcircle.fellenius import Slice, fellenius_fs

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    surf = CircularSurface(3.189, 16.909, 18.091)
    print(f"{'dovelas':>8} {'FS Bishop':>10} {'iteraciones':>12}")
    for n in (10, 20, 40, 80, 160, 320):
        r = BishopAnalysis(slope, n_slices=n).factor_of_safety(surf)
        print(f"{n:8d} {r['FS']:10.4f} {r.get('iterations', '?'):>12}")

    sl = slope.discretize_circle(surf, n_slices=160)
    fs_fell = fellenius_fs([Slice(s["weight_kN"], degrees(s["alpha_rad"]), s["cohesion"], degrees(s["phi_rad"]),
                                  s["width"] / cos(s["alpha_rad"]), s["u_kPa"]) for s in sl])
    fs_bishop = BishopAnalysis(slope, n_slices=160).factor_of_safety(surf)["FS"]
    print(f"\nFellenius (160 dovelas) = {fs_fell:.4f} ≤ Bishop = {fs_bishop:.4f}: {fs_fell <= fs_bishop} "
          f"(diferencia {fs_bishop / fs_fell - 1:+.1%})")


if __name__ == "__main__":
    main()
