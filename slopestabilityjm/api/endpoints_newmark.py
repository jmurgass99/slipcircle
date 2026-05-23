"""Endpoint Newmark sliding block (mejora SS-003 expuesta como REST)."""
from __future__ import annotations


def register_newmark_endpoint(app):
    from pydantic import BaseModel, Field

    from ..newmark import NewmarkInput, evaluate_newmark

    class NewmarkPayload(BaseModel):
        accelerogram_g: list[float] = Field(..., min_length=10)
        dt_s: float = Field(..., gt=0)
        a_y_g: float = Field(..., gt=0)

    @app.post("/api/slope/newmark")
    def newmark(payload: NewmarkPayload):
        r = evaluate_newmark(NewmarkInput(**payload.model_dump()))
        return {
            "permanent_displacement_m": r.permanent_displacement_m,
            "permanent_displacement_cm": r.permanent_displacement_m * 100,
            "n_pulses_exceeding": r.n_pulses_exceeding,
            "arias_intensity_m_s": r.arias_intensity_m_s,
            "severity": r.severity,
        }
