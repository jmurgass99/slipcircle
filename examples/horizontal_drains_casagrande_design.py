"""Drenes horizontales en un talud: número de drenes y tiempo de descenso (ecuación de Casagrande simplificada).

Talud de 12 m con permeabilidades de 10⁻⁷ a 10⁻⁵ m/s, drenes de 15 m cada 5 m, objetivo
de bajar el nivel 3 m. `design_drains` da el caudal por dren Q = k·H·L·f y el número de drenes
para evacuar el caudal objetivo; `time_to_drawdown_s` estima el tiempo de descenso a una
distancia dada (Theis simplificado). Son fórmulas de anteproyecto: el proyecto de drenaje se
hace con red de flujo o modelo numérico.

Ejecutar:  python examples/horizontal_drains_casagrande_design.py
"""
from slipcircle.drains_casagrande import DrainDesignInput, design_drains, time_to_drawdown_s


def main() -> None:
    print(f"{'k [m/s]':>8} {'Q/dren [l/min]':>15} {'n drenes':>9} {'Q total [m³/día]':>17} {'t descenso a 20 m':>18}")
    for k in (1e-7, 1e-6, 1e-5):
        r = design_drains(DrainDesignInput(k_m_s=k, slope_height_m=12.0, target_drawdown_m=3.0,
                                           drain_length_m=15.0, drain_spacing_m=5.0, factor=0.6, code="EU"))
        t_days = time_to_drawdown_s(k, 3.0, 20.0, porosity=0.3) / 86400.0
        print(f"{k:8.0e} {r.Q_per_drain_m3_s * 60000:15.3f} {r.n_drains_required:9d} {r.total_Q_m3_day:17.3f} {t_days:15.1f} días")
    r = design_drains(DrainDesignInput(k_m_s=1e-6, slope_height_m=12.0, target_drawdown_m=3.0,
                                       drain_length_m=25.0, drain_spacing_m=5.0))
    print(f"\nCon drenes de 25 m: {r.n_drains_required} drenes por tramo de 5 m, Q total {r.total_Q_m3_day:.3f} m³/día")


if __name__ == "__main__":
    main()
