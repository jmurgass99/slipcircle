"""Geotextil de refuerzo en la base de un terraplén sobre arcilla: incremento de FS como momento resistente en Bishop.

La ficha W00222 («Mejora del terreno», UPC, §6) da la cláusula: la fuerza T del geotextil
entra como momento resistente en el análisis de Bishop. Terraplén 1V:2H de 10 m sobre arcilla
no drenada (c_u = 25 kPa, γ = 18), círculo (15, 22, 24) y geotextil horizontal en y = 0 con
T_adm = 50 kN/m y τ_adherencia = 30 kPa. Para φ' = 0 el incremento tiene forma cerrada:
ΔFS = T·(c_y − y_geo) / Σ W·(x_m − c_x) = 50·22 / Σ W·(x_m − c_x). La longitud anclada
(24,6 m) da un arrancamiento de 1 476 kN/m, muy por encima de T_adm, que es la que gobierna.

Ejecutar:  python examples/geotextile_at_embankment_base_bishop_moment_w00222.py
"""
from slipcircle import BishopAnalysis, CircularSurface, GeosyntheticLayer, Slope, SoilLayer, reinforcement_moment

CX, CY, R, N = 15.0, 22.0, 24.0, 200


def main() -> None:
    slope = Slope(profile=[(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)],
                  layers=[SoilLayer("arcilla", cohesion=25.0, friction_angle=0.0, unit_weight=18.0)],
                  layer_depths=[40.0])
    surf = CircularSurface(CX, CY, R)
    geo = GeosyntheticLayer(y=0.0, x_left=0.0, x_right=30.0, T_allow=50.0, tau_bond=30.0)

    fs0 = BishopAnalysis(slope, n_slices=N).factor_of_safety(surf)["FS"]
    r1 = BishopAnalysis(slope, n_slices=N, reinforcements=[geo]).factor_of_safety(surf)
    slices = slope.discretize_circle(surf, N)
    driving = sum(s["weight_kN"] * (s["x_mid"] - CX) for s in slices)
    x_int = geo.intersect_circle(surf)
    T = geo.design_force(x_int, sigma_n=0.0)
    print(f"Intersección geotextil-círculo en x = {x_int:.3f} m; fuerza de diseño T = {T:.1f} kN/m (gobierna T_adm)")
    print(f"Momento del refuerzo (reinforcement_moment) = {reinforcement_moment([geo], surf, slope):.1f} kN·m/m")
    print(f"FS sin refuerzo = {fs0:.4f}")
    print(f"FS con geotextil = {r1['FS']:.4f}  (M_refuerzo en el resultado: {r1.get('M_reinforcement_kNm', '?')})")
    print(f"ΔFS cerrado = T·(c_y − y)/Σ W·(x_m − c_x) = {T * (CY - 0.0) / driving:.4f}; ΔFS calculado = {r1['FS'] - fs0:.4f}")


if __name__ == "__main__":
    main()
