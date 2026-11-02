"""Terraplén sobre arcilla blanda: rotura profunda por la base y efecto de una sobrecarga en coronación.

Terraplén de 6 m (1V:1,5H) de relleno granular (c' = 2 kPa, φ' = 32°, γ = 20) sobre arcilla
blanda no drenada (c_u = 20 kPa, φ = 0, γ = 17). La búsqueda por defecto obliga a que el
círculo salga por el pie; con `require_toe_exit=False` y `max_depth_below_toe=1.0` se
permiten círculos profundos que cortan la arcilla, que son los críticos en este caso. Después
se añade una sobrecarga uniforme de 20 kPa en coronación (`Slope.surcharge`).

La interfaz relleno/arcilla es horizontal (y = 0) bajo todo el terraplén, así que el talud se
define con `layer_reference="elevation"` y `layer_depths=[0, 20]` (base del relleno en y = 0,
base de la arcilla en y = -20). Con el valor por defecto (`"surface"`) los espesores se miden
desde la superficie local y la interfaz quedaría paralela al perfil: bajo el pie habría 6 m de
relleno en lugar de arcilla. El ejemplo lo comprueba.

Ejecutar:  python examples/embankment_on_soft_clay_deep_failure_and_surcharge.py
"""
from slipcircle import BishopAnalysis, CircularSurface, CriticalSurfaceSearch, Slope, SoilLayer

H = 6.0
PROFILE = [(0.0, 0.0), (1.5 * H, H), (1.5 * H + 20.0, H)]
LAYERS = [SoilLayer("relleno granular", cohesion=2.0, friction_angle=32.0, unit_weight=20.0),
          SoilLayer("arcilla blanda", cohesion=20.0, friction_angle=0.0, unit_weight=17.0)]
DEPTHS = [0.0, 20.0]   # profundidades bajo y = 0 (layer_reference="elevation")


def search(slope: Slope, **kw):
    res = CriticalSurfaceSearch(slope, n_cx=14, n_cy=12, n_R=12, n_slices=30, **kw).search()
    cs = res["critical_surface"]
    return res["critical_FS"], CircularSurface(cs["cx"], cs["cy"], cs["R"]), res["n_valid"]

def main() -> None:
    slope = Slope(profile=PROFILE, layers=LAYERS, layer_depths=DEPTHS, layer_reference="elevation")
    parallel = Slope(profile=PROFILE, layers=LAYERS, layer_depths=[H, 20.0])
    print(f"Estrato en (x = -3, y = -2), delante del pie: elevation → {slope.layer_at_depth(-2.0, x=-3.0).name}; "
          f"surface con [H, 20] → {parallel.layer_at_depth(-2.0, x=-3.0).name}")
    fs_toe, s_toe, n1 = search(slope)
    fs_deep, s_deep, n2 = search(slope, require_toe_exit=False, max_depth_below_toe=1.0)
    print(f"Solo círculos que salen por el pie:      FS = {fs_toe:.3f}  centro ({s_toe.cx:.1f}, {s_toe.cy:.1f}) R = {s_toe.R:.1f}  ({n1} válidos)")
    print(f"Círculos profundos permitidos (base):    FS = {fs_deep:.3f}  centro ({s_deep.cx:.1f}, {s_deep.cy:.1f}) R = {s_deep.R:.1f}  ({n2} válidos)")
    y_min = s_deep.cy - s_deep.R
    print(f"  cota mínima del círculo profundo: {y_min:+.2f} m (la arcilla empieza en 0 m)")
    cu, gamma_f = LAYERS[1].cohesion, LAYERS[0].unit_weight
    ref = 5.14 * cu / (gamma_f * H)
    print(f"  referencia a mano (hundimiento de la base, φ = 0): N_c·c_u/(γ·H) = 5,14·{cu:.0f}/({gamma_f:.0f}·{H:.0f}) = {ref:.2f}"
          f" → el círculo profundo da {fs_deep / ref:.2f} veces esa referencia")

    loaded = Slope(profile=PROFILE, layers=LAYERS, layer_depths=DEPTHS, layer_reference="elevation", surcharge=20.0)
    fs_q = BishopAnalysis(loaded, n_slices=40).factor_of_safety(s_deep)["FS"]
    fs_0 = BishopAnalysis(slope, n_slices=40).factor_of_safety(s_deep)["FS"]
    print(f"\nSobre el círculo profundo: FS sin sobrecarga = {fs_0:.3f}, con 20 kPa en coronación = {fs_q:.3f} ({fs_q / fs_0 - 1:+.1%})")


if __name__ == "__main__":
    main()
