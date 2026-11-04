# slipcircle

**Español** · [English](README.en.md)

Estabilidad de taludes por equilibrio límite: Fellenius, Bishop, Janbu, Spencer y Morgenstern-Price.

[![CI](https://github.com/jmurgass99/slipcircle/actions/workflows/ci.yml/badge.svg)](https://github.com/jmurgass99/slipcircle/actions/workflows/ci.yml)
[![Licencia: MIT](https://img.shields.io/badge/licencia-MIT-green)](LICENSE)

## Qué hace

- Factor de seguridad por el método de las dovelas sobre superficies circulares y poligonales, y búsqueda del círculo crítico.
- Presión de poros por freático, línea piezométrica o r_u, y sobrecarga en coronación.
- Pseudo-estático con los coeficientes de EN 1998-5 y desplazamiento permanente (Newmark, Bray-Travasarou).
- Refuerzo con geosintéticos y bulones, y drenes horizontales de anteproyecto.
- Los métodos iterativos devuelven `converged`: si es `False`, el FS es la mejor estimación y no un resultado válido.

## Instalación

```bash
pip install git+https://github.com/jmurgass99/slipcircle.git
```

Extras opcionales: `api`, `ui`. Dependencias base: `numpy` y `scipy`.

## Ejemplo

Talud de 10 m con pendiente 30° en una arena arcillosa (c = 10 kPa, φ = 20°, γ = 18 kN/m³):

```python
from math import tan, radians
from slipcircle import Slope, SoilLayer, CircularSurface, BishopAnalysis

dx = 10.0 / tan(radians(30.0))
suelo = SoilLayer(name="arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)
talud = Slope(profile=[(0.0, 0.0), (dx, 10.0), (dx + 30.0, 10.0)],
              layers=[suelo], layer_depths=[30.0])

r = BishopAnalysis(talud, n_slices=40).factor_of_safety(CircularSurface(cx=3.189, cy=16.909, R=18.091))
print(round(r["FS"], 3))   # 1.362
```

Más casos en [`examples/`](examples).

## Qué no hace

- Redes de flujo ni consolidación: la presión de poros se da, no se calcula.
- Roturas por cuña en 3D ni taludes rocosos (no incluye Hoek-Brown).
- Diseño EC7 con factores parciales: compara valores característicos con mínimos tabulados.
- Empuje sísmico en muros (Mononobe-Okabe): `mononobe` es una cuña pseudo-estática de talud indefinido.

El detalle está en [docs/LIMITACIONES.md](docs/LIMITACIONES.md).

## Pruebas

```bash
pip install -e ".[test]"
pytest -q
```

## Licencia

MIT. Ver [LICENSE](LICENSE). Autor: [Jesús Murgas](https://github.com/jmurgass99).
