r"""Fase J noche 2 (2026-09-16): comparacion contra patch paralelo independiente
(Auditoria paralela de contratos HTTP) encontro una inconsistencia real en `api_mononobe`.

Todo el resto de endpoints de este router que fallan por payload incompleto
devuelven 422 via `_wrap()` (que ya mapea KeyError/TypeError/ValueError desde
el commit 9bd9ac9). `api_mononobe` es el UNICO que interceptaba KeyError antes
de llegar a `_wrap()` y devolvia 400 con un formato de mensaje distinto -- una
frontera HTTP inconsistente, no un false-pass, pero si un contrato irregular
que un cliente no puede predecir sin leer caso por caso.
"""
from __future__ import annotations

from fastapi.testclient import TestClient

from slipcircle.api.app import app

client = TestClient(app)


def test_mononobe_payload_incompleto_devuelve_422_no_400() -> None:
    res = client.post("/api/slope/mononobe", json={})
    assert res.status_code == 422, (
        f"el resto del router usa 422 para payload incompleto via _wrap(); "
        f"mononobe debe ser consistente, no un segundo contrato: {res.text}"
    )
