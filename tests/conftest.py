"""Configuracion compartida de pytest.

- Asegura que el paquete ``slopestabilityjm`` sea importable cuando los
  tests se ejecutan desde el directorio raiz sin instalacion editable.
- Fixtures comunes (``simple_slope``, ``test_surface``) reutilizables.
"""
from __future__ import annotations
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from slopestabilityjm.soil_profile import SoilLayer  # noqa: E402
from slopestabilityjm.geometry import Slope, CircularSurface  # noqa: E402


@pytest.fixture
def simple_slope() -> Slope:
    """Talud 1V:2H, H=10 m, c'=20 kPa, phi'=20 deg, gamma=18 kN/m3."""
    layer = SoilLayer("Soil", cohesion=20.0, friction_angle=20.0, unit_weight=18.0)
    return Slope(
        profile=[(0, 0), (20, 10), (30, 10)],
        layers=[layer], layer_depths=[20.0],
    )


@pytest.fixture
def test_surface() -> CircularSurface:
    """Superficie circular tipica para el ``simple_slope``."""
    return CircularSurface(cx=15.0, cy=22.0, R=24.0)
