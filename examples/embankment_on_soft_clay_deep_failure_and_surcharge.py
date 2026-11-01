"""Terraplén sobre arcilla blanda: rotura profunda por la base y efecto de una sobrecarga en coronación.

Terraplén de 6 m (1V:1,5H) de relleno granular (c' = 2 kPa, φ' = 32°, γ = 20) sobre arcilla
blanda no drenada (c_u = 20 kPa, φ = 0, γ = 17). La búsqueda por defecto obliga a que el
círculo salga por el pie; con `require_toe_exit=False` y `max_depth_below_toe=1.0` se
permiten círculos profundos que cortan la arcilla, que son los críticos en este caso. Después
se añade una sobrecarga uniforme de 20 kPa en coronación (`Slope.surcharge`).

Ejecutar:  python examples/embankment_on_soft_clay_deep_failure_and_surcharge.py
"""
from slipcircle import BishopAnalysis, CircularSurface, CriticalSurfaceSearch, Slope, SoilLayer

H = 6.0
PROFILE = [(0.0, 0.0), (1.5 * H, H), (1.5 * H + 20.0, H)]
LAYERS = [SoilLayer("relleno granular", cohesion=2.0, friction_angle=32.0, unit_weight=20.0),
          SoilLayer("arcilla blanda", cohesion=20.0, friction_angle=0.0, unit_weight=17.0)]
DEPTHS = [H, 20.0]


def search(slope: Slope, **kw):
    res = CriticalSurfaceSearch(slope, n_cx=14, n_cy=12, n_R=12, n_slices=30, **kw).search()
    cs = res["critical_surface"]
    return res["critical_FS"], CircularSurface(cs["cx"], cs["cy"], cs["R"]), res["n_valid"]

def main() -> None:
    slope = Slope(profile=PROFILE, layers=LAYERS, layer_depths=DEPTHS)
    fs_toe, s_toe, n1 = search(slope)
    fs_deep, s_deep, n2 = search(slope, require_toe_exit=False, max_depth_below_toe=1.0)
    print(f"Solo círculos que salen por el pie:      FS = {fs_toe:.3f}  centro ({s_toe.cx:.1f}, {s_toe.cy:.1f}) R = {s_toe.R:.1f}  ({n1} válidos)")
    print(f"Círculos profundos permitidos (base):    FS = {fs_deep:.3f}  centro ({s_deep.cx:.1f}, {s_deep.cy:.1f}) R = {s_deep.R:.1f}  ({n2} válidos)")
    y_min = s_deep.cy - s_deep.R
    print(f"  cota mínima del círculo profundo: {y_min:+.2f} m (la arcilla empieza en 0 m)")

    loaded = Slope(profile=PROFILE, layers=LAYERS, layer_depths=DEPTHS, surcharge=20.0)
    fs_q = BishopAnalysis(loaded, n_slices=40).factor_of_safety(s_deep)["FS"]
    fs_0 = BishopAnalysis(slope, n_slices=40).factor_of_safety(s_deep)["FS"]
    print(f"\nSobre el círculo profundo: FS sin sobrecarga = {fs_0:.3f}, con 20 kPa en coronación = {fs_q:.3f} ({fs_q / fs_0 - 1:+.1%})")


if __name__ == "__main__":
    main()
