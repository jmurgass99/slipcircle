"""Búsqueda del círculo crítico: cuánto cambia el FS con la finura de la malla de centros y radios.

Talud homogéneo de Duncan y Wright (H = 10 m, β = 30°, c' = 10 kPa, φ' = 20°, γ = 18 kN/m³).
`CriticalSurfaceSearch` barre centros (n_cx × n_cy) y radios (n_R) con Bishop y devuelve el
mínimo. La referencia publicada para este talud está entre 1,30 y 1,40 según método y búsqueda.
El ejemplo imprime el FS mínimo de cada malla y lo compara con esa horquilla: el mínimo no tiene
por qué bajar al refinar, porque cada malla muestrea centros y radios distintos.

Ejecutar:  python examples/critical_surface_search_grid_refinement_duncan_wright.py
"""
import time
from math import radians, tan

from slipcircle import CriticalSurfaceSearch, Slope, SoilLayer

H, BETA = 10.0, radians(30.0)
GRIDS = [(8, 8, 6), (15, 12, 12), (24, 18, 14)]


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    print(f"{'malla (cx, cy, R)':<18} {'superficies':>11} {'válidas':>8} {'FS mín':>8} {'cx':>7} {'cy':>7} {'R':>7} {'s':>6}")
    fs_by_grid = []
    for ncx, ncy, nR in GRIDS:
        t0 = time.perf_counter()
        res = CriticalSurfaceSearch(slope, n_cx=ncx, n_cy=ncy, n_R=nR, n_slices=30).search()
        dt = time.perf_counter() - t0
        cs = res["critical_surface"]
        fs_by_grid.append(res["critical_FS"])
        print(f"{str((ncx, ncy, nR)):<18} {res['n_surfaces_tried']:11d} {res['n_valid']:8d} {res['critical_FS']:8.4f} "
              f"{cs['cx']:7.2f} {cs['cy']:7.2f} {cs['R']:7.2f} {dt:6.1f}")
    print("\nLas cinco mejores superficies de la última malla:")
    for r in res["all_results"][:5]:
        print(f"  FS = {r['FS']:.4f}  centro ({r['cx']:.2f}, {r['cy']:.2f})  R = {r['R']:.2f}")
    lo, hi = min(fs_by_grid), max(fs_by_grid)
    print(f"\nFS mínimo entre mallas: {lo:.3f} – {hi:.3f}; horquilla publicada 1,30 – 1,40 → "
          f"{'dentro' if lo >= 1.30 else f'{(1.30 - lo) / 1.30:.1%} por debajo del límite inferior'}"
          f"{'' if lo >= 1.30 or hi < 1.30 else ' en la malla más desfavorable; el resto dentro'}")
    print(f"acceptable_EC7 (FS ≥ 1,5 con parámetros característicos): {res['acceptable_EC7']}")


if __name__ == "__main__":
    main()
