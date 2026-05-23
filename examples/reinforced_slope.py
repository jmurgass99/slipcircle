"""
Reinforced version of the Duncan & Wright homogeneous slope.

Same geometry and soil as ``duncan_wright_homogeneous.py`` but with a single
geosynthetic (geogrid) layer at mid-height.  The reinforcement moment is wired
directly into each LEM, so the user only needs to pass ``reinforcements=[...]``
to the analysis constructor — no manual ΔMr post-processing required.

Run with:

    python -m examples.reinforced_slope
"""
from math import tan, radians

from slopestabilityjm import (
    Slope, SoilLayer, CircularSurface,
    BishopAnalysis, JanbuAnalysis, SpencerAnalysis,
    MorgensternPriceAnalysis,
    GeosyntheticLayer,
)


def build_slope() -> Slope:
    H = 10.0
    beta = radians(30.0)
    dx = H / tan(beta)
    profile = [(0.0, 0.0), (dx, H), (dx + 30.0, H)]
    soil = SoilLayer(
        name="Homogeneous clayey sand",
        cohesion=10.0, friction_angle=20.0, unit_weight=18.0,
    )
    return Slope(profile=profile, layers=[soil], layer_depths=[30.0])


def main() -> None:
    slope = build_slope()
    # Critical circle from the unreinforced search (see duncan_wright example).
    surf = CircularSurface(cx=3.189, cy=16.909, R=18.091)

    # Geogrid at mid-height (y = 5 m) spanning the slope face plus 5 m beyond
    # the crest as anchorage length.
    H = 10.0
    dx_face = H / tan(radians(30.0))
    geogrid = GeosyntheticLayer(
        y=5.0,
        x_left=0.5 * dx_face,         # under the slope face
        x_right=dx_face + 5.0,        # 5 m behind the crest
        T_allow=80.0,                 # kN/m design tensile capacity
        tau_bond=35.0,                # kPa
    )

    print("=== FS comparison: unreinforced vs reinforced ===")
    for cls in [BishopAnalysis, JanbuAnalysis, SpencerAnalysis,
                MorgensternPriceAnalysis]:
        r0 = cls(slope, n_slices=40).factor_of_safety(surf)
        r1 = cls(slope, n_slices=40,
                 reinforcements=[geogrid]).factor_of_safety(surf)
        delta = r1["FS"] - r0["FS"]
        print(f"  {r0['method']:30s}  "
              f"FS_no_reinf={r0['FS']:.4f}  "
              f"FS_reinf={r1['FS']:.4f}  "
              f"dFS=+{delta:.4f}")


if __name__ == "__main__":
    main()
