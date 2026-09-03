"""Tests de integración y golden para análisis pseudoestático Mononobe.

Norma / Referencia: Mononobe (1926) / Eurocódigo 8 Parte 5 Anexo E (Fichas W00678/W00704).
"""
import math
import pytest
from fastapi.testclient import TestClient

from slopestabilityjm.api.app import app
from slopestabilityjm.mononobe import MononobeInput, evaluate_mononobe


@pytest.fixture
def client():
    return TestClient(app)


def test_mononobe_golden_hand_calc():
    """Golden analítico exacto derivado a mano (W00678 / W00704).

    Talud α = 30°, H = 5.0 m, γ = 18.0 kN/m³, c = 10.0 kPa, φ = 30°:
      - FS estático:
          c/(γ·H·sinα·cosα) + tanφ/tanα = 10/(18·5·0.5·0.866025) + 1.0 = 0.25660 + 1.0 = 1.2566
      - Sismo kh = 0.2, kv = 0.1:
          ψ = atan(0.2 / 0.9) = 12.5288°
    """
    inp = MononobeInput(
        slope_height_m=5.0,
        slope_angle_deg=30.0,
        gamma_kN_m3=18.0,
        c_kPa=10.0,
        phi_deg=30.0,
        kh=0.2,
        kv=0.1,
    )
    res = evaluate_mononobe(inp, FS_min=1.1)

    assert res.FS_static == pytest.approx(1.2566, abs=1e-3)
    assert res.psi_deg == pytest.approx(12.5288, abs=1e-3)
    assert res.FS_seismic < res.FS_static


def test_mononobe_api_endpoint(client):
    """Verifica endpoint HTTP /api/slope/mononobe."""
    payload = {
        "slope_height_m": 5.0,
        "slope_angle_deg": 30.0,
        "gamma_kN_m3": 18.0,
        "c_kPa": 10.0,
        "phi_deg": 30.0,
        "kh": 0.2,
        "kv": 0.1,
        "FS_min": 1.1,
    }
    resp = client.post("/api/slope/mononobe", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["FS_static"] == pytest.approx(1.2566, abs=1e-3)
    assert data["psi_deg"] == pytest.approx(12.5288, abs=1e-3)
    assert "FS_seismic" in data
    assert data["is_safe"] == (data["FS_seismic"] >= 1.1)
