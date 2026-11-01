"""Factor de seguridad mínimo por norma y veredicto: EC7 (DA con γ_R), USACE/FHWA, NSR-10, AS 4678.

El mismo FS calculado (Bishop sobre el talud de Duncan y Wright) se compara con el mínimo
que exige cada código en situación estática persistente y sísmica. En Eurocódigo 7 el
mínimo es 1,0 porque los factores parciales ya van en los parámetros (aquí se usan valores
característicos, así que la comparación con 1,0 no es un diseño EC7 completo); en US, CO y
AU se comparan valores característicos con 1,5 / 1,05–1,1.

Ejecutar:  python examples/fs_required_by_code_and_verdict.py
"""
from math import radians, tan

from slipcircle import BishopAnalysis, CircularSurface, PseudoStaticAnalysis, Slope, SoilLayer
from slipcircle.fs_by_code import FS_MIN_TABLE, fs_minimum, is_safe

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    surf = CircularSurface(3.189, 16.909, 18.091)
    fs_static = BishopAnalysis(slope, n_slices=40).factor_of_safety(surf)["FS"]
    fs_seismic = PseudoStaticAnalysis(slope=slope, kh=0.10, kv=0.05, n_slices=40).factor_of_safety(surf)["FS"]
    print(f"FS estático (Bishop) = {fs_static:.3f};  FS pseudo-estático (k_h = 0,10, k_v = 0,05) = {fs_seismic:.3f}\n")
    print(f"{'código':<6} {'escenario':<18} {'FS mín':>7} {'FS calc':>8}  veredicto")
    for (code, scenario), fs_min in FS_MIN_TABLE.items():
        fs = fs_static if scenario == "static_persistent" else fs_seismic
        print(f"{code:<6} {scenario:<18} {fs_min:7.2f} {fs:8.3f}  {'cumple' if is_safe(fs, code, scenario) else 'NO cumple'}")
    try:
        fs_minimum("EU", "accidental")
    except ValueError as exc:
        print(f"\nEscenario no tabulado → ValueError: {exc}")


if __name__ == "__main__":
    main()
