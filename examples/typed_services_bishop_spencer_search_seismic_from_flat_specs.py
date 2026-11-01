"""Servicios tipados: el mismo talud como especificación plana serializable (JSON) y cuatro análisis.

`SlopeSpec`, `SoilLayerSpec` y `CircleSpec` son dataclasses inmutables que se guardan o se
envían por HTTP sin esfuerzo; los servicios `bishop`, `spencer`, `slip_search` y
`seismic_slope` las convierten y devuelven resultados tipados (con el desplazamiento de
Newmark si se da `ag_peak`). Una geometría inválida lanza las excepciones del paquete.

Ejecutar:  python examples/typed_services_bishop_spencer_search_seismic_from_flat_specs.py
"""
from dataclasses import asdict
from math import radians, tan

from slipcircle.exceptions import GeometryError
from slipcircle.services import bishop, seismic_slope, slip_search, spencer
from slipcircle.services._common import CircleSpec, SlopeSpec, SoilLayerSpec

H, BETA = 10.0, radians(30.0)


def main() -> None:
    dx = H / tan(BETA)
    slope = SlopeSpec(profile=((0.0, 0.0), (dx, H), (dx + 30.0, H)),
                      layers=(SoilLayerSpec("arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0, depth_bottom=30.0),))
    circle = CircleSpec(cx=3.189, cy=16.909, R=18.091)

    print("Bishop:  ", asdict(bishop.analyze(bishop.BishopInput(slope=slope, circle=circle, n_slices=40))))
    print("Spencer: ", asdict(spencer.analyze(spencer.SpencerInput(slope=slope, circle=circle, n_slices=40))))
    print("Búsqueda:", asdict(slip_search.analyze(slip_search.SlipSearchInput(slope=slope, n_cx=12, n_cy=10, n_R=10))))
    print("Sísmico: ", asdict(seismic_slope.analyze(seismic_slope.SeismicSlopeInput(slope=slope, circle=circle, kh=0.10, kv=0.05, ag_peak=0.25))))
    try:
        bishop.analyze(bishop.BishopInput(slope=slope, circle=CircleSpec(cx=3.0, cy=16.0, R=-1.0)))
    except GeometryError as exc:
        print(f"\nRadio negativo → GeometryError: {exc}")


if __name__ == "__main__":
    main()
