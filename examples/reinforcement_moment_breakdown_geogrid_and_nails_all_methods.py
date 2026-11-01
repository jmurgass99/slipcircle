"""Momento resistente de geomalla y bulones combinados, y cómo lo recogen los cuatro métodos.

Talud de Duncan y Wright con una geomalla a media altura (T_adm = 80 kN/m) y dos bulones
de 10 m. `reinforcement_moment` suma T por su brazo respecto al centro del círculo (vertical
para la geomalla horizontal, perpendicular a la barra para el bulón); Bishop, Janbu, Spencer
y Morgenstern-Price aceptan la misma lista `reinforcements` y el FS sube de forma coherente
entre métodos.

Ejecutar:  python examples/reinforcement_moment_breakdown_geogrid_and_nails_all_methods.py
"""
from math import radians, tan

from slipcircle import (
    BishopAnalysis,
    CircularSurface,
    GeosyntheticLayer,
    JanbuAnalysis,
    MorgensternPriceAnalysis,
    Slope,
    SoilLayer,
    SoilNail,
    SpencerAnalysis,
    reinforcement_moment,
)

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    surf = CircularSurface(3.189, 16.909, 18.091)
    geogrid = GeosyntheticLayer(y=5.0, x_left=0.5 * dx, x_right=dx + 5.0, T_allow=80.0, tau_bond=35.0)
    nails = [SoilNail(x_head=y / tan(BETA), y_head=y, inclination_deg=15.0, length=10.0, T_allow=150.0) for y in (3.0, 7.0)]

    for label, elems in (("geomalla", [geogrid]), ("2 bulones", nails), ("geomalla + 2 bulones", [geogrid] + nails)):
        print(f"M_refuerzo {label:<22} = {reinforcement_moment(elems, surf, slope):8.1f} kN·m/m")

    all_elems = [geogrid] + nails
    print(f"\n{'método':<22} {'FS sin':>8} {'FS con':>8} {'ΔFS':>7}")
    for cls in (BishopAnalysis, JanbuAnalysis, SpencerAnalysis, MorgensternPriceAnalysis):
        r0 = cls(slope, n_slices=40).factor_of_safety(surf)
        r1 = cls(slope, n_slices=40, reinforcements=all_elems).factor_of_safety(surf)
        print(f"{r0['method']:<22} {r0['FS']:8.4f} {r1['FS']:8.4f} {r1['FS'] - r0['FS']:+7.4f}")


if __name__ == "__main__":
    main()
