"""
SlopeStabilityJm — Slope Stability Analysis Library
====================================================
Methods: Bishop, Janbu, Spencer, Morgenstern-Price, Pseudo-static seismic.

Usage
-----
>>> from slopestabilityjm import Slope, SoilLayer, BishopAnalysis
>>> slope = Slope(...)
>>> result = BishopAnalysis(slope).factor_of_safety()
"""
from .geometry import Slope, CircularSurface
from .soil_profile import SoilLayer
from .slices import slices_from_circle, slices_from_polyline
from .bishop import BishopAnalysis
from .janbu import JanbuAnalysis
from .janbu_generalized import JanbuGeneralizedAnalysis
from .spencer import SpencerAnalysis
from .morgenstern_price import MorgensternPriceAnalysis, half_sine, constant_f, trapezoidal
from .seismic import PseudoStaticAnalysis
from .search import CriticalSurfaceSearch
from .water import HorizontalWaterTable, PiezometricLine, RuModel
from .reinforcement import GeosyntheticLayer, SoilNail, reinforcement_moment
from .exceptions import (
    SlopeError, GeotechError, ConvergenceError, CodeError,
    GeometryError, ConfigurationError, DependencyError,
)

__all__ = [
    "Slope", "CircularSurface", "SoilLayer",
    "slices_from_circle", "slices_from_polyline",
    "BishopAnalysis", "JanbuAnalysis", "JanbuGeneralizedAnalysis",
    "SpencerAnalysis", "MorgensternPriceAnalysis",
    "half_sine", "constant_f", "trapezoidal",
    "PseudoStaticAnalysis", "CriticalSurfaceSearch",
    "HorizontalWaterTable", "PiezometricLine", "RuModel",
    "GeosyntheticLayer", "SoilNail", "reinforcement_moment",
    # excepciones tipadas
    "SlopeError", "GeotechError", "ConvergenceError", "CodeError",
    "GeometryError", "ConfigurationError", "DependencyError",
]
__version__ = "1.2.0"
