"""Dibujar el talud, el círculo crítico y las dovelas a PNG (extra `ui`, matplotlib sin ventana).

Talud de Duncan y Wright con su círculo crítico y 40 dovelas; `plot_result` superpone perfil,
círculo etiquetado con el FS y las dovelas. Se guarda en examples/output/. Requiere
`pip install "slipcircle[ui]"`.

Ejecutar:  python examples/plot_critical_circle_and_slices_png.py
"""
import os
from math import radians, tan

import matplotlib

matplotlib.use("Agg")

from slipcircle import BishopAnalysis, CircularSurface, Slope, SoilLayer, slices_from_circle  # noqa: E402
from slipcircle.visualization import plot_result  # noqa: E402

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = Slope(profile=[(0.0, 0.0), (dx, H), (dx + 30.0, H)],
                  layers=[SoilLayer("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)],
                  layer_depths=[30.0])
    surf = CircularSurface(3.189, 16.909, 18.091)
    res = BishopAnalysis(slope, n_slices=40).factor_of_safety(surf)
    ax = plot_result(slope, surf, res, slices=slices_from_circle(slope, surf, 40), title="Duncan & Wright, Bishop")
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "duncan_wright_critical_circle.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    ax.figure.savefig(out, dpi=110, bbox_inches="tight")
    print(f"FS = {res['FS']}; figura guardada en {os.path.relpath(out, os.getcwd())} ({os.path.getsize(out) // 1024} kB)")


if __name__ == "__main__":
    main()
