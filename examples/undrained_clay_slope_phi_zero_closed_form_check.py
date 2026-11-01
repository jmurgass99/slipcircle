"""Talud en arcilla no drenada (φ = 0): Bishop y Fellenius frente a la forma cerrada del círculo.

Talud 1V:2H de 10 m en arcilla con c_u = 25 kPa y γ = 18 kN/m³, círculo de centro
(15, 22) y radio 24 m. Para φ' = 0 los dos métodos de dovelas colapsan a la expresión
exacta FS = c·R·L_arco / Σ W·(x_m − c_x): el momento resistente es el de la cohesión
sobre todo el arco. Es el anclaje analítico de `tests/test_golden_slip_methods.py`
(el método de Bishop sobre círculos lo cita la ficha W00254 del máster).

Ejecutar:  python examples/undrained_clay_slope_phi_zero_closed_form_check.py
"""
from math import asin, cos, degrees

from slipcircle import BishopAnalysis, CircularSurface, Slope, SoilLayer
from slipcircle.fellenius import Slice, fellenius_fs

C_U, CX, CY, R, N = 25.0, 15.0, 22.0, 24.0, 200


def main() -> None:
    slope = Slope(profile=[(0.0, 0.0), (20.0, 10.0), (30.0, 10.0)],
                  layers=[SoilLayer("arcilla", cohesion=C_U, friction_angle=0.0, unit_weight=18.0)],
                  layer_depths=[40.0])
    surf = CircularSurface(CX, CY, R)
    slices = slope.discretize_circle(surf, N)

    x_lo = slices[0]["x_mid"] - slices[0]["width"] / 2
    x_hi = slices[-1]["x_mid"] + slices[-1]["width"] / 2
    L_arc = R * abs(asin((x_hi - CX) / R) - asin((x_lo - CX) / R))
    driving = sum(s["weight_kN"] * (s["x_mid"] - CX) for s in slices)
    fs_closed = C_U * R * L_arc / driving

    fs_bishop = BishopAnalysis(slope, n_slices=N).factor_of_safety(surf)["FS"]
    fs_fell = fellenius_fs([Slice(s["weight_kN"], degrees(s["alpha_rad"]), s["cohesion"], degrees(s["phi_rad"]),
                                  s["width"] / cos(s["alpha_rad"]), s["u_kPa"]) for s in slices])
    print(f"Arco entre x = {x_lo:.2f} y {x_hi:.2f} m: L_arco = {L_arc:.3f} m; Σ W·(x_m − c_x) = {driving:.1f} kN·m/m")
    print(f"Forma cerrada c·R·L_arco/Σ W·(x_m − c_x) = {fs_closed:.4f}")
    print(f"Bishop simplificado ({N} dovelas)          = {fs_bishop:.4f}")
    print(f"Fellenius (Σ c·l, l = b/cos α)             = {fs_fell:.4f}")
    print("\nEn φ' = 0 los tres coinciden: cualquier diferencia sería un error en el tratamiento de mα o de la base de la dovela.")


if __name__ == "__main__":
    main()
