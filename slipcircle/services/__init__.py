"""Service layer de SlipCircle.

Funciones puras tipadas (dataclass Input + Result) que envuelven los
metodos LEM (Bishop, Janbu, Spencer), busqueda de superficie critica y
analisis pseudo-estatico sismico. Reutilizables desde GUI, dashboard,
REST API y job queue sin acoplamiento a la clase ``Slope`` original.
"""
from __future__ import annotations
from . import bishop, janbu, spencer, slip_search, seismic_slope

__all__ = ["bishop", "janbu", "spencer", "slip_search", "seismic_slope"]
