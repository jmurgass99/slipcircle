"""Bray y Travasarou (2007): desplazamiento sísmico mediano y probabilidades de excedencia.

Dos casos anclados a la ecuación publicada (Bray 2007, cap. 14, ec. 14.5, idéntica a
Bray & Travasarou 2007): bloque casi rígido (T_s = 0: k_y = 0,10, Sa = 0,40 g, M = 7,0 →
15,4 cm con el intercepto −0,22) y el caso de validación Austrian Dam en Loma Prieta
(Tabla 14.1: k_y = 0,14, T_s = 0,33 s, Sa(1,5·T_s) = 0,94 g, M = 6,93 → 34,2 cm con
−1,10 + 1,5·T_s). Hasta agosto de 2026 el módulo mezclaba los dos interceptos y daba
82 cm para Austrian Dam. Las probabilidades P(D > 5 cm) y P(D > 30 cm) salen de la
lognormal con σ = 0,66.

Ejecutar:  python examples/bray_travasarou_2007_rigid_block_and_austrian_dam.py
"""
from slipcircle.bray_travasarou import BrayTravasarouInput, evaluate_bray_travasarou

CASES = [
    ("bloque rígido (T_s = 0)", BrayTravasarouInput(k_y=0.10, Sa_1_5Ts_g=0.40, Mw=7.0, Ts_seconds=0.0), "15,4 cm"),
    ("Austrian Dam, Loma Prieta 1989", BrayTravasarouInput(k_y=0.14, Sa_1_5Ts_g=0.94, Mw=6.93, Ts_seconds=0.33), "34,2 cm (observado 20–70 cm)"),
    ("mismo talud con k_y = 0,25", BrayTravasarouInput(k_y=0.25, Sa_1_5Ts_g=0.94, Mw=6.93, Ts_seconds=0.33), "—"),
]


def main() -> None:
    print(f"{'caso':<32} {'D mediano':>10} {'P(D>5cm)':>9} {'P(D>30cm)':>10}  referencia")
    for name, inp, ref in CASES:
        r = evaluate_bray_travasarou(inp)
        print(f"{name:<32} {r.median_displacement_cm:8.1f} cm {r.P_D_gt_5cm:9.3f} {r.P_D_gt_30cm:10.3f}  {ref}")
    print(f"\nσ_lnD = {evaluate_bray_travasarou(CASES[0][1]).sigma_lnD}; nivel de validación declarado por el módulo: "
          f"{evaluate_bray_travasarou(CASES[0][1]).validation_level}")


if __name__ == "__main__":
    main()
