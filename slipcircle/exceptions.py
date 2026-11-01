"""Jerarquia de excepciones tipadas de SlipCircle.

Permite al caller (REST API, GUI, scripts) diferenciar la
naturaleza del fallo y reaccionar apropiadamente.

Todas las clases heredan tambien de ``ValueError`` para conservar la
compatibilidad hacia atras con codigo que captura ``ValueError`` generico
(la libreria historicamente lanzaba ``ValueError`` para entradas invalidas).

Ejemplo
-------
>>> from slipcircle.exceptions import GeotechError, ConvergenceError
>>> try:
...     analysis.factor_of_safety(surface)
... except ConvergenceError as exc:
...     print(exc.iteration, exc.residual)
... except GeotechError:
...     # parametros de suelo invalidos
...     ...
"""
from __future__ import annotations


class SlopeError(ValueError):
    """Base de todas las excepciones del paquete.

    Hereda de ``ValueError`` para mantener compatibilidad con codigo
    cliente existente que captura ``ValueError`` generico.
    """


class GeotechError(SlopeError):
    """Parametros geotecnicos invalidos o inconsistentes.

    Ejemplos: phi negativo o >=90, cohesion negativa, peso unitario nulo,
    perfil con espesores incoherentes, nivel freatico sobre la cresta.
    """


class ConvergenceError(SlopeError):
    """El solver iterativo (Bishop, Janbu, Spencer, MP, Newmark) no convergio.

    Atributos
    ---------
    iteration: int
        Ultima iteracion alcanzada antes de abortar.
    residual: float
        Residuo numerico al fallar (|FS - FS_prev| o |sumF|+|sumM|).
    """

    def __init__(self, message: str, *, iteration: int = -1,
                 residual: float = float("nan")) -> None:
        super().__init__(message)
        self.iteration = iteration
        self.residual = residual


class CodeError(SlopeError):
    """Caso fuera del rango de aplicabilidad de la norma seleccionada.

    Ejemplos: EC8-5 con kh > coeficiente maximo permitido, EC7 con
    geometria fuera del Annex B, Mononobe-Okabe con phi-theta-beta <= 0.
    """


class GeometryError(SlopeError):
    """Geometria de talud o superficie de falla invalida.

    Ejemplos: superficie circular que no intersecta el terreno, polilinea
    no monotona en x, perfil con menos de 2 puntos, momento motor no
    positivo (superficie fisicamente inviable).
    """


class ConfigurationError(SlopeError):
    """Parametro de configuracion invalido o incompatible."""


class DependencyError(SlopeError):
    """Falta un paquete opcional (fastapi, hypothesis, matplotlib, ...).

    El mensaje debe indicar el extra a instalar, e.g.
    ``pip install slipcircle[api]``.
    """


__all__ = [
    "SlopeError",
    "GeotechError",
    "ConvergenceError",
    "CodeError",
    "GeometryError",
    "ConfigurationError",
    "DependencyError",
]
