"""
Classic Duncan & Wright homogeneous slope example.

  Height H = 10 m, inclination β = 30° (1V : sqrt(3)H ≈ 1V : 1.732H)
  c' = 10 kPa,  φ' = 20°,  γ = 18 kN/m³,  dry

Expected FS (literature): ~1.30-1.40 depending on the method and the
critical-circle search.  Bishop tends to give ~1.34; rigorous methods
(Spencer, Morgenstern-Price) are within ±2%.

Run with:

    python -m examples.duncan_wright_homogeneous
"""
from math import tan, radians

from slipcircle import (
    Slope, SoilLayer, CircularSurface,
    BishopAnalysis, JanbuAnalysis, SpencerAnalysis,
    MorgensternPriceAnalysis, CriticalSurfaceSearch,
    PseudoStaticAnalysis,
)


def build_slope() -> Slope:
    H = 10.0
    beta = radians(30.0)
    dx = H / tan(beta)
    # toe at (0, 0); slope face up to crest at (dx, H); crest plateau to right.
    profile = [(0.0, 0.0), (dx, H), (dx + 30.0, H)]
    soil = SoilLayer(
        name="Homogeneous clayey sand",
        cohesion=10.0, friction_angle=20.0, unit_weight=18.0,
    )
    return Slope(profile=profile, layers=[soil], layer_depths=[30.0])


def main() -> None:
    slope = build_slope()

    print("=== Critical surface search (Bishop) ===")
    search = CriticalSurfaceSearch(slope, n_cx=15, n_cy=12, n_R=12, n_slices=30)
    res = search.search()
    cs = res["critical_surface"]
    print(f"  Critical FS = {res['critical_FS']}")
    print(f"  Critical circle: cx={cs['cx']}  cy={cs['cy']}  R={cs['R']}")

    surf = CircularSurface(cs["cx"], cs["cy"], cs["R"])

    print("\n=== Compare methods on the critical circle ===")
    for analyzer in [
        BishopAnalysis(slope, n_slices=40),
        JanbuAnalysis(slope, n_slices=40),
        SpencerAnalysis(slope, n_slices=40),
        MorgensternPriceAnalysis(slope, n_slices=40),
    ]:
        r = analyzer.factor_of_safety(surf)
        print(f"  {r['method']:30s}  FS = {r.get('FS')}")

    print("\n=== Pseudo-static (kh = 0.10) ===")
    pse = PseudoStaticAnalysis(slope=slope, kh=0.10, kv=0.05, n_slices=40)
    print("  ", pse.factor_of_safety(surf))


if __name__ == "__main__":
    main()
