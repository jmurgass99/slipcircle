"""Pseudo-estático según EN 1998-5 §4.1.3.3: k_h = α·S/(2·r), k_v = 0,5·k_h, y desplazamiento estimado.

Talud de Duncan y Wright sobre su círculo crítico para a_g = 0,10 / 0,15 / 0,25 g, factor de
suelo S = 1,0 y 1,2, y r = 1,0 (desplazamientos pequeños) o 1,5 (hasta 20 cm admitidos).
`PseudoStaticAnalysis.from_ec8` construye los coeficientes; después `newmark_displacement`
da el desplazamiento permanente por la fórmula empírica de Ambraseys y Menu (1988) con
k_y = k_h y k_max = a_g, y lo compara con el límite de 30 cm.

Ejecutar:  python examples/pseudostatic_ec8_5_kh_from_ag_soil_factor_and_newmark_estimate.py
"""
from math import radians, tan

from slipcircle import CircularSurface, PseudoStaticAnalysis, Slope, SoilLayer

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    surf = CircularSurface(3.189, 16.909, 18.091)
    print(f"{'a_g [g]':>7} {'S':>4} {'r':>4} {'k_h':>6} {'k_v':>6} {'FS':>7}  Newmark (Ambraseys-Menu)")
    for ag in (0.10, 0.15, 0.25):
        for S, r in ((1.0, 1.0), (1.2, 1.0), (1.2, 1.5)):
            ps = PseudoStaticAnalysis.from_ec8(slope, ag=ag, S=S, r=r, n_slices=40)
            fs = ps.factor_of_safety(surf)["FS"]
            nm = ps.newmark_displacement(FS_yield=1.0, ag_peak=ag * S)
            d = nm.get("displacement_m", 0.0)
            print(f"{ag:7.2f} {S:4.1f} {r:4.1f} {ps.kh:6.3f} {ps.kv:6.3f} {fs:7.3f}  D = {100 * d:5.1f} cm, "
                  f"admisible (< 30 cm): {nm.get('acceptable', nm.get('note', '?'))}")
    print("\nAmbraseys-Menu usa k_y = k_h del análisis, no la aceleración de fluencia real del talud (FS = 1): "
          "es una estimación de cribado, no un cálculo de bloque deslizante.")


if __name__ == "__main__":
    main()
