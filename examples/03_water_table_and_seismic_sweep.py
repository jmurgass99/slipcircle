"""Influencia del nivel freático y del coeficiente sísmico en el factor de seguridad.

Talud de 12 m, 1V:1,5H, dos estratos (relleno arenoso sobre arcilla). Se busca
la superficie crítica en seco, se repite con el nivel freático a 4 m sobre el
pie y después se hace un barrido pseudo-estático de k_h (EC8-5) sobre la
superficie crítica del caso húmedo.

Ejecutar:  python examples/03_water_table_and_seismic_sweep.py
"""
from slipcircle import (
    BishopAnalysis,
    CircularSurface,
    CriticalSurfaceSearch,
    PseudoStaticAnalysis,
    Slope,
    SoilLayer,
)

H = 12.0
PROFILE = [(0.0, 0.0), (1.5 * H, H), (1.5 * H + 25.0, H)]
LAYERS = [
    SoilLayer(name="relleno arenoso", cohesion=5.0, friction_angle=30.0, unit_weight=19.0),
    SoilLayer(name="arcilla", cohesion=25.0, friction_angle=18.0, unit_weight=18.0),
]
DEPTHS = [6.0, 30.0]


def critical(slope: Slope):
    res = CriticalSurfaceSearch(slope, n_cx=12, n_cy=10, n_R=10, n_slices=30).search()
    cs = res["critical_surface"]
    return res["critical_FS"], CircularSurface(cs["cx"], cs["cy"], cs["R"])


def main() -> None:
    dry = Slope(profile=PROFILE, layers=LAYERS, layer_depths=DEPTHS)
    wet = Slope(profile=PROFILE, layers=LAYERS, layer_depths=DEPTHS, water_table_y=4.0)

    fs_dry, surf_dry = critical(dry)
    fs_wet, surf_wet = critical(wet)
    print(f"FS crítico en seco:            {fs_dry:.3f}")
    print(f"FS crítico con freático a 4 m: {fs_wet:.3f}  (Δ = {100 * (fs_wet / fs_dry - 1):+.1f} %)\n")

    print("Barrido pseudo-estático sobre la superficie crítica húmeda (Bishop modificado):")
    print(f"  {'k_h':>5} {'FS':>7}")
    for kh in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25):
        if kh == 0.0:
            fs = BishopAnalysis(wet, n_slices=40).factor_of_safety(surf_wet)["FS"]
        else:
            fs = PseudoStaticAnalysis(slope=wet, kh=kh, kv=0.0, n_slices=40).factor_of_safety(surf_wet)["FS"]
        flag = "" if fs >= 1.0 else "  ← inestable"
        print(f"  {kh:5.2f} {fs:7.3f}{flag}")


if __name__ == "__main__":
    main()
