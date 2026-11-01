"""API REST (extra `api`) sin servidor: Bishop, círculo crítico, pseudo-estático, cuña sísmica y Newmark por HTTP.

Los endpoints reciben el talud como JSON (perfil, capas con `depth_bottom`, freático,
sobrecarga) y el círculo como {cx, cy, R}. Una entrada incompleta devuelve 422. La aplicación
se sirve con `slipcircle api` o `uvicorn slipcircle.api.app:app`.

Ejecutar:  python examples/rest_api_bishop_critical_slip_seismic_newmark_testclient.py
"""
import json
import warnings
from math import radians, tan

warnings.filterwarnings("ignore")
from fastapi.testclient import TestClient  # noqa: E402

from slipcircle.api.app import app  # noqa: E402

H, BETA = 10.0, radians(30.0)


def main() -> None:
    c = TestClient(app)
    dx = H / tan(BETA)
    slope = {"profile": [[0.0, 0.0], [dx, H], [dx + 30.0, H]],
             "layers": [{"name": "arena arcillosa", "cohesion": 10.0, "friction_angle": 20.0, "unit_weight": 18.0, "depth_bottom": 30.0}]}
    circle = {"cx": 3.189, "cy": 16.909, "R": 18.091}

    print("GET /api/slope/health →", c.get("/api/slope/health").json())
    for path, body in (("/api/slope/bishop", {"slope": slope, "circle": circle, "n_slices": 40}),
                       ("/api/slope/critical_slip", {"slope": slope, "n_cx": 10, "n_cy": 8, "n_R": 8}),
                       ("/api/slope/seismic", {"slope": slope, "circle": circle, "kh": 0.10, "kv": 0.05, "ag_peak": 0.25}),
                       ("/api/slope/mononobe", {"slope_height_m": 5.0, "slope_angle_deg": 30.0, "gamma_kN_m3": 18.0,
                                                "c_kPa": 10.0, "phi_deg": 30.0, "kh": 0.2, "kv": -0.1})):
        r = c.post(path, json=body)
        print(f"\nPOST {path} → {r.status_code}")
        print("  ", json.dumps(r.json(), ensure_ascii=False)[:400])

    acc = [0.5] * 500 + [0.0] * 1000
    r = c.post("/api/slope/newmark", json={"accelerogram_g": acc, "dt_s": 0.001, "a_y_g": 0.2})
    print(f"\nPOST /api/slope/newmark (pulso 0,5 g × 0,5 s, k_y = 0,2) → {r.status_code}", json.dumps(r.json(), ensure_ascii=False) if r.status_code == 200 else r.text[:200])

    r = c.post("/api/slope/bishop", json={"slope": slope})
    print(f"\nSin círculo → {r.status_code}")


if __name__ == "__main__":
    main()
