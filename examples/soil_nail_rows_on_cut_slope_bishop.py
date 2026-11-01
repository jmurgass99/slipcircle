"""Filas de bulones (soil nails) en un desmonte: fuerza de diseño por anclaje o por barra y FS con Bishop.

Desmonte de Duncan y Wright con una, dos y tres filas de bulones de 10 m a 15° de inclinación,
cabezas en la cara del talud a 2,5, 5,0 y 7,5 m, T_adm = 150 kN, adherencia 80 kPa, taladro
de 100 mm. Para cada bulón se calcula la longitud que queda tras la superficie de rotura y
la fuerza que gobierna (arrancamiento π·d·τ·L o capacidad de la barra); Bishop suma su
momento respecto al centro del círculo.

Ejecutar:  python examples/soil_nail_rows_on_cut_slope_bishop.py
"""
from math import radians, tan

from slipcircle import BishopAnalysis, CircularSurface, Slope, SoilLayer, SoilNail

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    surf = CircularSurface(3.189, 16.909, 18.091)
    nails = [SoilNail(x_head=y / tan(BETA), y_head=y, inclination_deg=15.0, length=10.0, T_allow=150.0,
                      bond_stress=80.0, diameter=0.10) for y in (2.5, 5.0, 7.5)]
    print("Bulón   cabeza (x, y)   intersección (x, y)   L tras rotura   T arrancamiento   T diseño")
    for n in nails:
        hit = n.intersect_circle(surf)
        if hit is None:
            print(f"  ({n.x_head:5.2f}, {n.y_head:4.1f})  no corta la superficie")
            continue
        xi, yi, L = hit
        T_po = 3.141592653589793 * n.diameter * n.bond_stress * L
        print(f"  ({n.x_head:5.2f}, {n.y_head:4.1f})   ({xi:5.2f}, {yi:5.2f})   {L:9.2f} m   {T_po:12.1f} kN   {n.design_force(L):8.1f} kN")

    fs0 = BishopAnalysis(slope, n_slices=40).factor_of_safety(surf)["FS"]
    print(f"\nFS sin bulones = {fs0:.4f}")
    for k in (1, 2, 3):
        r = BishopAnalysis(slope, n_slices=40, reinforcements=nails[:k]).factor_of_safety(surf)
        print(f"FS con {k} fila(s) = {r['FS']:.4f}  (ΔFS = {r['FS'] - fs0:+.4f})")


if __name__ == "__main__":
    main()
