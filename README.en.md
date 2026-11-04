# slipcircle

[Español](README.md) · **English**

Slope stability by limit equilibrium: Fellenius, Bishop, Janbu, Spencer and Morgenstern-Price.

[![CI](https://github.com/jmurgass99/slipcircle/actions/workflows/ci.yml/badge.svg)](https://github.com/jmurgass99/slipcircle/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## What it does

- Factor of safety by the method of slices on circular and polygonal surfaces, and critical-circle search.
- Pore pressure from a water table, a piezometric line or r_u, and surcharge at the crest.
- Pseudo-static analysis with the EN 1998-5 coefficients and permanent displacement (Newmark, Bray-Travasarou).
- Reinforcement with geosynthetics and soil nails, and preliminary horizontal drains.
- The iterative methods return `converged`: if it is `False`, the FS is the best estimate and not a valid result.

## Install

```bash
pip install git+https://github.com/jmurgass99/slipcircle.git
```

Optional extras: `api`, `ui`. Base dependencies: `numpy` and `scipy`.

## Example

A 10 m slope at 30° in a clayey sand (c = 10 kPa, φ = 20°, γ = 18 kN/m³):

```python
from math import tan, radians
from slipcircle import Slope, SoilLayer, CircularSurface, BishopAnalysis

dx = 10.0 / tan(radians(30.0))
soil = SoilLayer(name="clayey sand", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)
slope = Slope(profile=[(0.0, 0.0), (dx, 10.0), (dx + 30.0, 10.0)],
              layers=[soil], layer_depths=[30.0])

r = BishopAnalysis(slope, n_slices=40).factor_of_safety(CircularSurface(cx=3.189, cy=16.909, R=18.091))
print(round(r["FS"], 3))   # 1.362
```

More cases in [`examples/`](examples).

## What it does not do

- Flow nets or consolidation: pore pressure is an input, not a result.
- 3D wedge failures or rock slopes (no Hoek-Brown).
- EC7 design with partial factors: it compares characteristic values with tabulated minimums.
- Seismic thrust on walls (Mononobe-Okabe): `mononobe` is an infinite-slope pseudo-static wedge.

Details are in [docs/LIMITACIONES.md](docs/LIMITACIONES.md) (in Spanish).

## Tests

```bash
pip install -e ".[test]"
pytest -q
```

## License

MIT. See [LICENSE](LICENSE). Author: [Jesús Murgas](https://github.com/jmurgass99).
