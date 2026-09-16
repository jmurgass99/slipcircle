"""Test de contrato HTTP 422 para errores de validación de payload en SlopeStabilityJm."""
import pytest
from fastapi.testclient import TestClient
from slopestabilityjm.api.app import app


@pytest.fixture
def client():
    return TestClient(app)


def test_bishop_incomplete_payload_returns_422(client):
    """Un payload incompleto debe devolver HTTP 422, nunca HTTP 500."""
    r1 = client.post("/api/slope/bishop", json={})
    assert r1.status_code == 422, f"Esperaba 422 pero obtuvo {r1.status_code}: {r1.text}"

    r2 = client.post("/api/slope/bishop", json={"slope": {}})
    assert r2.status_code == 422, f"Esperaba 422 pero obtuvo {r2.status_code}: {r2.text}"


def test_janbu_incomplete_payload_returns_422(client):
    """Janbu con payload incompleto devuelve HTTP 422."""
    r = client.post("/api/slope/janbu", json={"slope": {"layers": []}})
    assert r.status_code == 422, f"Esperaba 422 pero obtuvo {r.status_code}: {r.text}"
