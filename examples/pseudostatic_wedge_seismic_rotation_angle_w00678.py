"""Cuña pseudo-estática de talud indefinido: ángulo de rotación sísmica ψ = atan(k_h/(1 ± k_v)) de las fichas W00678/W00704.

`mononobe.evaluate_mononobe` NO es el K_AE de Mononobe-Okabe para muros: es una cuña de talud
indefinido en la que la gravedad aparente gira ψ. Las fichas «Seismic design of foundations»
(Prof. A. Drougkas, UPC) imprimen ψ = 10,3° para k_h = 0,2 con +k_v = 0,1 y 12,5° con −k_v;
el módulo define ψ = atan(k_h/(1 − k_v)), así que el caso «+k_v» de la ficha se obtiene con
k_v = −0,1. Talud α = 30°, H = 5 m, γ = 18, c' = 10 kPa, φ' = 30°. El FS estático coincide
con la forma cerrada del talud indefinido.

Ejecutar:  python examples/pseudostatic_wedge_seismic_rotation_angle_w00678.py
"""
from math import cos, radians, sin, tan

from slipcircle.mononobe import MononobeInput, evaluate_mononobe

BASE = dict(slope_height_m=5.0, slope_angle_deg=30.0, gamma_kN_m3=18.0, c_kPa=10.0, phi_deg=30.0)


def main() -> None:
    a = radians(30.0)
    fs_closed = (10.0 / cos(a) + 18.0 * 5.0 * cos(a) * tan(radians(30.0))) / (18.0 * 5.0 * sin(a))
    print(f"FS estático cerrado (talud indefinido) = {fs_closed:.4f}\n")
    print(f"{'k_h':>5} {'k_v':>6} {'ψ [°]':>6} {'FS estático':>11} {'FS sísmico':>10}  seguro (FS ≥ 1,1)")
    for kh, kv, label in ((0.0, 0.0, ""), (0.2, -0.1, "ficha «+k_v»: 10,3°"), (0.2, 0.1, "ficha «−k_v»: 12,5°"),
                          (0.1, 0.0, ""), (0.3, 0.0, "")):
        r = evaluate_mononobe(MononobeInput(**BASE, kh=kh, kv=kv))
        print(f"{kh:5.2f} {kv:+6.2f} {r.psi_deg:6.2f} {r.FS_static:11.4f} {r.FS_seismic:10.4f}  {r.is_safe}  {label}")


if __name__ == "__main__":
    main()
