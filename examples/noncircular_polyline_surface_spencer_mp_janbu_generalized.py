"""Superficie de rotura no circular (polilínea) en un talud estratificado: cuatro métodos sobre la misma polilínea.

Talud de 12 m, 1V:1,5H, relleno granular (c' = 5 kPa, φ' = 30°) de 6 m sobre arcilla
(c' = 25 kPa, φ' = 18°). La superficie pasa por la arcilla y sale por la coronación: una
polilínea de cinco vértices que ningún círculo reproduce. Janbu generalizado itera las fuerzas
entre dovelas; Janbu simplificado da la cota rápida con f₀; Spencer y Morgenstern-Price toman el
equilibrio de momentos respecto al punto medio de la cuerda y devuelven `converged` (FS de
fuerzas = FS de momentos) junto con los dos FS. Si `converged` fuese False, el FS devuelto no
debe usarse.

Ejecutar:  python examples/noncircular_polyline_surface_spencer_mp_janbu_generalized.py
"""
from slipcircle import (
    JanbuAnalysis,
    JanbuGeneralizedAnalysis,
    MorgensternPriceAnalysis,
    Slope,
    SoilLayer,
    SpencerAnalysis,
)

H = 12.0
PROFILE = [(0.0, 0.0), (1.5 * H, H), (1.5 * H + 25.0, H)]
LAYERS = [SoilLayer("relleno granular", cohesion=5.0, friction_angle=30.0, unit_weight=19.0),
          SoilLayer("arcilla", cohesion=25.0, friction_angle=18.0, unit_weight=18.0)]
POLYLINE = [(0.0, 0.0), (6.0, -1.0), (14.0, 0.5), (22.0, 6.5), (28.0, 12.0)]


def main() -> None:
    slope = Slope(profile=PROFILE, layers=LAYERS, layer_depths=[6.0, 30.0])
    print("Superficie polilínea:", POLYLINE, "\n")
    r = JanbuGeneralizedAnalysis(slope, surface=POLYLINE, n_slices=40).compute_fs()
    print(f"Janbu generalizado:    FS = {r['FS']:.4f}  iteraciones = {r.get('iterations', '?')}  convergió = {r['converged']}")
    r = JanbuAnalysis(slope, n_slices=40).factor_of_safety(POLYLINE)
    print(f"Janbu simplificado:    FS = {r['FS']:.4f}  f₀ = {r['correction_f0']}")
    r = SpencerAnalysis(slope, n_slices=40).factor_of_safety(POLYLINE)
    print(f"Spencer (polilínea):   FS = {r['FS']:.4f}  (fuerzas {r['FS_force_eq']:.4f}, momentos {r['FS_moment_eq']:.4f})  "
          f"θ = {r['theta_deg']:.1f}°  convergió = {r['converged']}")
    r = MorgensternPriceAnalysis(slope, n_slices=40).factor_of_safety(POLYLINE)
    print(f"M-P (polilínea):       FS = {r['FS']:.4f}  (fuerzas {r['FS_force_eq']:.4f}, momentos {r['FS_moment_eq']:.4f})  "
          f"λ = {r['lambda']:.2f}  convergió = {r['converged']}")
    print("\nJanbu generalizado queda por debajo de Spencer y M-P (no impone la dirección de las fuerzas "
          "entre dovelas y usa la línea de empuje por defecto); Janbu simplificado con f₀ es la cota rápida.")


if __name__ == "__main__":
    main()
