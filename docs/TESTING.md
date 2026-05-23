# Testing — SlopeStabilityJm

## Layout

- `tests/test_slope_stability.py` — unit tests historicos (Bishop, Janbu, Spencer, MP, seismic, search, refuerzos).
- `tests/test_exceptions.py` — jerarquia tipada `SlopeError → GeotechError | ConvergenceError | CodeError | GeometryError | ConfigurationError | DependencyError` y modo `strict=True` en Bishop.
- `tests/test_services.py` — service layer (`slopestabilityjm.services.*`), puente FoundationDesignJm y REST API (FastAPI TestClient).
- `tests/test_property_invariants.py` — invariantes fisicos verificados con [hypothesis](https://hypothesis.readthedocs.io/) (FS positivo, monotonia en c'/phi'/kh, FS sismico <= FS estatico, Bishop ~ Janbu).
- `tests/conftest.py` — fixtures `simple_slope` y `test_surface`.

## Markers (registrados en `pyproject.toml`)

| Marker            | Significado                                    |
| ----------------- | ---------------------------------------------- |
| `unit`            | unitarios rapidos (< 100 ms)                   |
| `integration`     | requieren stack completo o multiples modulos   |
| `slow`            | tests con duracion alta (> 5 s)                |
| `safety_critical` | validados contra valores publicados (tol <= 2%) |
| `property`        | property-based con hypothesis                  |
| `api`             | requieren FastAPI (`pip install -e .[api]`)    |

## Ejecucion

```powershell
# Suite completa
python -m pytest tests/ -q

# Solo property-based
python -m pytest tests/ -m property -q

# Con cobertura (aplica el gate fail_under=75% del pyproject)
.\scripts\run_coverage.ps1
```

## Gate de cobertura

`pyproject.toml` declara `[tool.coverage.report] fail_under = 75`. El script
`run_coverage.ps1` ejecuta `pytest --cov=slopestabilityjm`; si la cobertura
total cae por debajo del 75% el comando termina con exit code != 0.

Modulos omitidos del calculo (declarados en `[tool.coverage.run] omit`):

- `slopestabilityjm/_brand.py` — constantes
- `slopestabilityjm/visualization.py` — matplotlib (UI-only)
- `slopestabilityjm/__main__.py` — CLI

## Extras opcionales

```powershell
pip install -e .[test]   # pytest + coverage + hypothesis
pip install -e .[api]    # fastapi + uvicorn + httpx
pip install -e .[ui]     # matplotlib + streamlit
```
