"""Tres modelos de presión de poros sobre la misma superficie: freático horizontal, línea piezométrica y r_u.

Talud de 10 m (1V:2H) en suelo c' = 10 kPa, φ' = 25°, γ = 19 kN/m³, con su círculo crítico en
seco. `slices_from_circle` admite un modelo de agua con método u(x, y): freático horizontal en
la cota del pie (y = 0), línea piezométrica que baja hasta el pie y r_u = 0,25 constante (que necesita la altura
de cada dovela: `u_from_slice`). El FS se evalúa con el método Ordinario sobre esas dovelas y,
para el freático horizontal, también con Bishop a través de `Slope.water_table_y`.

Ejecutar:  python examples/pore_pressure_models_horizontal_piezometric_ru.py
"""
from math import cos, degrees

from slipcircle import (
    BishopAnalysis,
    CircularSurface,
    CriticalSurfaceSearch,
    HorizontalWaterTable,
    PiezometricLine,
    RuModel,
    Slope,
    SoilLayer,
    slices_from_circle,
)
from slipcircle.fellenius import Slice, fellenius_fs

PROFILE = [(0.0, 0.0), (20.0, 10.0), (40.0, 10.0)]
SOIL = SoilLayer("limo arenoso", cohesion=10.0, friction_angle=25.0, unit_weight=19.0)


def fs_fellenius(slices, u_override=None):
    out = []
    for s in slices:
        u = u_override(s) if u_override else s["u_kPa"]
        out.append(Slice(s["weight_kN"], degrees(s["alpha_rad"]), s["cohesion"], degrees(s["phi_rad"]),
                         s["width"] / cos(s["alpha_rad"]), u))
    return fellenius_fs(out)


def main() -> None:
    dry = Slope(profile=PROFILE, layers=[SOIL], layer_depths=[30.0])
    cs = CriticalSurfaceSearch(dry, n_cx=12, n_cy=10, n_R=10, n_slices=30).search()["critical_surface"]
    surf = CircularSurface(cs["cx"], cs["cy"], cs["R"])
    print(f"Círculo crítico en seco: centro ({cs['cx']}, {cs['cy']}), R = {cs['R']}\n")

    models = {
        "seco": None,
        "freático horizontal y_w = 0 m (pie)": HorizontalWaterTable(y_w=0.0),
        "línea piezométrica (0, 0) → (20, 6) → (40, 8)": PiezometricLine([(0.0, 0.0), (20.0, 6.0), (40.0, 8.0)]),
    }
    for name, model in models.items():
        sl = slices_from_circle(dry, surf, n_slices=60, water_model=model)
        u_max = max(s["u_kPa"] for s in sl)
        print(f"{name:<46} u_máx = {u_max:6.1f} kPa  FS Fellenius = {fs_fellenius(sl):.4f}")
    ru = RuModel(ru=0.25)
    sl = slices_from_circle(dry, surf, n_slices=60)
    print(f"{'r_u = 0,25 (u = r_u·γ·h por dovela)':<46} u_máx = {max(ru.u_from_slice(s) for s in sl):6.1f} kPa  "
          f"FS Fellenius = {fs_fellenius(sl, ru.u_from_slice):.4f}")

    wet = Slope(profile=PROFILE, layers=[SOIL], layer_depths=[30.0], water_table_y=0.0)
    print(f"\nBishop con Slope.water_table_y = 0 m: FS = {BishopAnalysis(wet, n_slices=60).factor_of_safety(surf)['FS']:.4f} "
          f"(seco: {BishopAnalysis(dry, n_slices=60).factor_of_safety(surf)['FS']:.4f})")


if __name__ == "__main__":
    main()
