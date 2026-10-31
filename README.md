# slipcircle

Estabilidad de taludes por equilibrio límite: Fellenius, Bishop simplificado, Janbu (simplificado y generalizado), Spencer y Morgenstern-Price, con búsqueda automática de la superficie crítica, nivel freático (freático horizontal, línea piezométrica, r_u), análisis pseudo-estático (EC8-5), desplazamiento de Newmark y Bray-Travasarou, y refuerzo con geosintéticos y anclajes.

[![CI](https://github.com/jmurgass99/slipcircle/actions/workflows/ci.yml/badge.svg)](https://github.com/jmurgass99/slipcircle/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

Paquete independiente. Dependencias: `numpy` y `scipy`.

## Instalación

```bash
pip install git+https://github.com/jmurgass99/slipcircle.git
```

Extras: `ui` (matplotlib y streamlit para visualizar el talud y las dovelas), `api` (REST con FastAPI), `test`.

```bash
git clone https://github.com/jmurgass99/slipcircle.git
cd slipcircle
pip install -e ".[test,api,ui]"
pytest -q
```

## Inicio rápido

Talud homogéneo clásico de Duncan y Wright: H = 10 m, β = 30°, c' = 10 kPa, φ' = 20°, γ = 18 kN/m³, seco. Con la malla de búsqueda de abajo el paquete da FS ≈ 1,25 (Bishop); la referencia publicada para este talud está entre 1,3 y 1,4 según el método y la finura de la búsqueda, así que afina `n_cx`, `n_cy` y `n_R` antes de comparar.

```python
from math import tan, radians
from slipcircle import Slope, SoilLayer, CircularSurface, BishopAnalysis, CriticalSurfaceSearch

H, beta = 10.0, radians(30.0)
dx = H / tan(beta)
profile = [(0.0, 0.0), (dx, H), (dx + 30.0, H)]          # pie en el origen, cara, corona
soil = SoilLayer(name="arena arcillosa", cohesion=10.0, friction_angle=20.0, unit_weight=18.0)
slope = Slope(profile=profile, layers=[soil], layer_depths=[30.0])

# Búsqueda de la superficie circular crítica con Bishop
res = CriticalSurfaceSearch(slope, n_cx=15, n_cy=12, n_R=12, n_slices=30).search()
print(res["critical_FS"], res["critical_surface"])   # ≈ 1.25 con esta malla

# Comparar métodos sobre esa misma superficie
cs = res["critical_surface"]
surf = CircularSurface(cs["cx"], cs["cy"], cs["R"])
print(BishopAnalysis(slope, n_slices=40).factor_of_safety(surf)["FS"])   # ≈ 1.25
```

Pseudo-estático con coeficientes sísmicos:

```python
from slipcircle import PseudoStaticAnalysis
PseudoStaticAnalysis(slope=slope, kh=0.10, kv=0.05, n_slices=40).factor_of_safety(surf)
```

El script `examples/duncan_wright_homogeneous.py` ejecuta este caso completo con los cuatro métodos. Hay también una interfaz de línea de comandos (`slipcircle --help`).

## Módulos

| Módulo | Contenido |
|---|---|
| `geometry`, `soil_profile`, `slices`, `water` | Perfil del talud, estratos, dovelas, presión de poros |
| `fellenius`, `bishop`, `janbu`, `janbu_generalized`, `spencer`, `morgenstern_price` | Métodos de equilibrio límite (funciones f(x) `half_sine`, `constant_f`, `trapezoidal`) |
| `search` | Barrido de centros y radios y refinamiento |
| `seismic`, `newmark`, `bray_travasarou`, `mononobe` | Pseudo-estático, desplazamiento permanente, empuje sísmico |
| `reinforcement`, `reinforced_slope`, `drains_casagrande` | Geosintéticos, anclajes, drenes |
| `fs_by_code` | Factores de seguridad exigidos por norma (EC7, otras) |
| `duncan_wright_homogeneous` | Ábacos de Duncan y Wright para talud homogéneo |
| `visualization` | Dibujo del talud, superficie y dovelas (extra `ui`) |
| `integration/foundationcalc_bridge` | Construcción del talud a partir de estratos de un diseño de cimentación (import diferido; opcional) |
| `api/` | REST API con FastAPI (extra `api`) |

## Validación

Los tests marcados `safety_critical` comparan con valores publicados con tolerancia ≤ 2 % (Duncan y Wright, ábacos y ejemplos de referencia); los tests `property` usan hypothesis para invariantes (por ejemplo, FS decrece al subir k_h). Consulta `docs/TESTING.md` y `docs/QC_CHECKLIST.md`.

## Licencia

MIT. Consulta [LICENSE](LICENSE).

---

Parte del toolkit estructural de código abierto **Loadpath**, de [Jesús Murgas](https://github.com/jmurgass99).
