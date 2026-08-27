"""Tests de la capa de servicios (Mejora A)."""
from __future__ import annotations
import pytest

from slopestabilityjm.services._common import (
    SoilLayerSpec, SlopeSpec, CircleSpec, build_slope, build_circle,
)
from slopestabilityjm.services import bishop, janbu, spencer, slip_search, seismic_slope
from slopestabilityjm.exceptions import GeometryError, GeotechError, CodeError


def _spec(c: float = 20.0, phi: float = 20.0,
          gamma: float = 18.0, water: float | None = None) -> SlopeSpec:
    return SlopeSpec(
        profile=((0.0, 0.0), (20.0, 10.0), (30.0, 10.0)),
        layers=(SoilLayerSpec("Soil", c, phi, gamma, depth_bottom=20.0),),
        water_table_y=water,
    )


def _circ() -> CircleSpec:
    return CircleSpec(cx=15.0, cy=22.0, R=24.0)


class TestCommon:

    def test_build_slope_ok(self):
        s = build_slope(_spec())
        assert s.height == pytest.approx(10.0)

    def test_build_slope_no_layers_raises(self):
        bad = SlopeSpec(profile=((0, 0), (1, 1)), layers=())
        with pytest.raises(GeotechError):
            build_slope(bad)

    def test_build_slope_short_profile_raises(self):
        bad = SlopeSpec(profile=((0, 0),),
                        layers=(SoilLayerSpec("X", 10, 20, 18, 5.0),))
        with pytest.raises(GeometryError):
            build_slope(bad)

    def test_build_circle_negative_R(self):
        with pytest.raises(GeometryError):
            build_circle(CircleSpec(0, 0, -1.0))


class TestBishopService:

    def test_analyze_returns_typed_result(self):
        r = bishop.analyze(bishop.BishopInput(slope=_spec(), circle=_circ(), n_slices=25))
        assert isinstance(r, bishop.BishopResult)
        assert r.FS > 0
        assert r.method == "Bishop Simplified"

    def test_invalid_surface_raises(self):
        r_in = bishop.BishopInput(slope=_spec(),
                                  circle=CircleSpec(cx=-100, cy=200, R=0.5),
                                  n_slices=10)
        with pytest.raises(GeometryError):
            bishop.analyze(r_in)


class TestJanbuService:
    def test_analyze(self):
        r = janbu.analyze(janbu.JanbuInput(slope=_spec(), circle=_circ()))
        assert r.FS > 0
        assert 1.0 <= r.correction_f0 <= 1.2


class TestSpencerService:
    def test_analyze(self):
        r = spencer.analyze(spencer.SpencerInput(slope=_spec(), circle=_circ(), n_slices=25))
        assert r.FS > 0


class TestSlipSearchService:
    def test_search(self):
        r = slip_search.analyze(slip_search.SlipSearchInput(
            slope=_spec(), n_cx=4, n_cy=4, n_R=4, n_slices=15))
        assert r.critical_FS > 0


class TestSeismicSlopeService:
    def test_seismic_reduces_fs(self):
        r_static = bishop.analyze(bishop.BishopInput(slope=_spec(), circle=_circ()))
        r_seis = seismic_slope.analyze(seismic_slope.SeismicSlopeInput(
            slope=_spec(), circle=_circ(), kh=0.15))
        assert r_seis.FS < r_static.FS

    def test_negative_kh_raises(self):
        with pytest.raises(CodeError):
            seismic_slope.analyze(seismic_slope.SeismicSlopeInput(
                slope=_spec(), circle=_circ(), kh=-0.1))

    def test_with_newmark(self):
        r = seismic_slope.analyze(seismic_slope.SeismicSlopeInput(
            slope=_spec(), circle=_circ(), kh=0.10, ag_peak=0.25))
        assert r.displacement_m is not None
        assert r.displacement_m >= 0


class TestIntegrationBridge:

    def test_layers_from_flat_objects(self):
        from slopestabilityjm.integration.foundationdesignjm_bridge import (
            layers_from_foundation_spec, _FlatLayer,
        )
        flats = [_FlatLayer(thickness=3.0, gamma=18.0, phi=25.0, c=5.0),
                 _FlatLayer(thickness=5.0, gamma=19.0, phi=0.0, qu=120.0)]
        layers, depths = layers_from_foundation_spec(flats)
        assert len(layers) == 2
        assert depths == [3.0, 8.0]
        # Capa 2: phi=0 + qu=120 -> c = qu/2 = 60
        assert layers[1].cohesion == pytest.approx(60.0)

    def test_slope_from_foundation_spec(self):
        from slopestabilityjm.integration.foundationdesignjm_bridge import (
            slope_from_foundation_spec, _FlatLayer,
        )
        flats = [_FlatLayer(thickness=20.0, gamma=18.0, phi=20.0, c=20.0)]
        slope = slope_from_foundation_spec(
            flats, profile=[(0, 0), (20, 10), (30, 10)])
        assert slope.height == pytest.approx(10.0)
        # Sanity: el slope construido produce FS valido con Bishop
        from slopestabilityjm.bishop import BishopAnalysis
        from slopestabilityjm.geometry import CircularSurface
        r = BishopAnalysis(slope, n_slices=20).factor_of_safety(
            CircularSurface(15, 22, 24))
        assert r["FS"] > 0

    def test_empty_layers_raises(self):
        from slopestabilityjm.integration.foundationdesignjm_bridge import (
            layers_from_foundation_spec,
        )
        with pytest.raises(GeotechError):
            layers_from_foundation_spec([])

    def test_negative_thickness_raises(self):
        from slopestabilityjm.integration.foundationdesignjm_bridge import (
            layers_from_foundation_spec, _FlatLayer,
        )
        with pytest.raises(GeotechError):
            layers_from_foundation_spec([_FlatLayer(thickness=-1.0, gamma=18.0)])


@pytest.mark.api
class TestRestAPI:
    """Tests del FastAPI app (requieren extra [api])."""

    def setup_method(self):
        pytest.importorskip("fastapi")
        from fastapi.testclient import TestClient
        from slopestabilityjm.api.app import app
        self.client = TestClient(app)

    def test_root(self):
        r = self.client.get("/")
        assert r.status_code == 200
        body = r.json()
        assert "endpoints" in body
        assert any("bishop" in e for e in body["endpoints"])

    def test_health(self):
        r = self.client.get("/api/slope/health")
        assert r.status_code == 200
        assert r.json()["status"] == "ok"

    def _payload(self) -> dict:
        return {
            "slope": {
                "profile": [[0, 0], [20, 10], [30, 10]],
                "layers": [{
                    "name": "Soil", "cohesion": 20.0,
                    "friction_angle": 20.0, "unit_weight": 18.0,
                    "depth_bottom": 20.0,
                }],
            },
            "circle": {"cx": 15.0, "cy": 22.0, "R": 24.0},
        }

    def test_bishop_endpoint(self):
        r = self.client.post("/api/slope/bishop", json=self._payload())
        assert r.status_code == 200, r.text
        body = r.json()
        assert body["FS"] > 0
        assert body["method"] == "Bishop Simplified"

    def test_janbu_endpoint(self):
        r = self.client.post("/api/slope/janbu", json=self._payload())
        assert r.status_code == 200
        assert r.json()["FS"] > 0

    def test_seismic_endpoint(self):
        p = self._payload()
        p["kh"] = 0.15
        r = self.client.post("/api/slope/seismic", json=p)
        assert r.status_code == 200
        assert r.json()["FS"] > 0

    def test_critical_slip_endpoint(self):
        p = {"slope": self._payload()["slope"], "n_cx": 4, "n_cy": 4, "n_R": 4, "n_slices": 15}
        r = self.client.post("/api/slope/critical_slip", json=p)
        assert r.status_code == 200
        assert r.json()["critical_FS"] > 0

    def test_bad_geometry_returns_422(self):
        p = self._payload()
        p["circle"] = {"cx": -100, "cy": 200, "R": 0.5}
        r = self.client.post("/api/slope/bishop", json=p)
        assert r.status_code == 422

    def test_negative_kh_returns_422(self):
        p = self._payload()
        p["kh"] = -0.1
        r = self.client.post("/api/slope/seismic", json=p)
        assert r.status_code == 422
