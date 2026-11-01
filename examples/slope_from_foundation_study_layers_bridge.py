"""Construir el talud desde los estratos de un estudio geotécnico de cimentación (puente a foundationcalc).

`slope_from_foundation_spec` acepta cualquier objeto con `thickness`, `gamma`, `phi`, `c` y `qu`
(duck typing: los `SoilLayerSpec` de foundationcalc o un dataclass propio) y crea el `Slope` con
los espesores acumulados; para una arcilla con q_u y sin c toma c_u = q_u/2. Sobre el talud
construido se busca el círculo crítico con Bishop.

Ejecutar:  python examples/slope_from_foundation_study_layers_bridge.py
"""
from dataclasses import dataclass

from slipcircle import CriticalSurfaceSearch
from slipcircle.exceptions import GeotechError
from slipcircle.integration.foundationcalc_bridge import layers_from_foundation_spec, slope_from_foundation_spec


@dataclass
class StudyLayer:
    thickness: float
    gamma: float
    phi: float = 0.0
    c: float = 0.0
    qu: float = 0.0


STUDY = [StudyLayer(3.0, 19.0, phi=30.0, c=2.0),          # relleno
         StudyLayer(5.0, 17.5, phi=0.0, qu=60.0),          # arcilla blanda: c_u = q_u/2 = 30 kPa
         StudyLayer(12.0, 20.0, phi=34.0, c=0.0)]          # arena densa


def main() -> None:
    layers, depths = layers_from_foundation_spec(STUDY, name_prefix="EST")
    for L, d in zip(layers, depths):
        print(f"  {L.name:<6} c = {L.cohesion:5.1f} kPa  φ = {L.friction_angle:4.1f}°  γ = {L.unit_weight:4.1f}  base a {d:4.1f} m")
    slope = slope_from_foundation_spec(STUDY, profile=[(0.0, 0.0), (12.0, 8.0), (40.0, 8.0)], water_table_y=2.0)
    res = CriticalSurfaceSearch(slope, n_cx=12, n_cy=10, n_R=10, n_slices=30, require_toe_exit=False,
                                max_depth_below_toe=1.0).search()
    cs = res["critical_surface"]
    print(f"\nTalud 8 m (1V:1,5H) con freático a 2 m: FS crítico = {res['critical_FS']:.3f}, "
          f"círculo ({cs['cx']}, {cs['cy']}) R = {cs['R']}, cota mínima {cs['cy'] - cs['R']:+.2f} m")
    try:
        layers_from_foundation_spec([StudyLayer(0.0, 19.0)])
    except GeotechError as exc:
        print(f"\nCapa de espesor 0 → GeotechError: {exc}")


if __name__ == "__main__":
    main()
