"""Puente desde FoundationDesignJm hacia SlipCircle.

Convierte una secuencia de ``SoilLayerSpec`` (formato plano de
FoundationDesignJm.services.bearing) en la lista de ``SoilLayer`` +
``layer_depths`` que entiende ``slipcircle.geometry.Slope``.

El objetivo es habilitar el flujo:

    FoundationDesignJm (perfil estratigrafico investigado para una
    cimentacion) -> SlipCircle (estabilidad global del talud
    sobre el que se asienta la cimentacion).

Import de FoundationDesignJm es lazy: si no esta instalado se acepta
cualquier objeto con los atributos requeridos (``thickness``, ``gamma``,
``phi``, ``c``, ``qu``) -- duck typing.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence, Any

from ..soil_profile import SoilLayer
from ..geometry import Slope
from ..exceptions import DependencyError, GeotechError


@dataclass(frozen=True)
class _FlatLayer:
    """Forma minima que debe tener una capa para ser puenteada."""
    thickness: float
    gamma: float
    phi: float = 0.0
    c: float = 0.0
    qu: float = 0.0


def layers_from_foundation_spec(
    fd_layers: Sequence[Any],
    *,
    name_prefix: str = "FD",
) -> tuple[list[SoilLayer], list[float]]:
    """Convierte capas FoundationDesignJm a ``(SoilLayer[], depth_bottom[])``.

    Parameters
    ----------
    fd_layers : Sequence
        Iterable de ``foundationcalc.services.bearing.SoilLayerSpec``
        (o cualquier objeto con atributos compatibles via duck typing).
    name_prefix : str
        Prefijo para nombrar las capas resultantes (``FD_1``, ``FD_2``...).

    Returns
    -------
    (layers, depths_bottom) listos para ``Slope(layers=..., layer_depths=...)``.
    """
    if not fd_layers:
        raise GeotechError("Sequence de capas vacia.")

    layers: list[SoilLayer] = []
    depths: list[float] = []
    acc = 0.0
    for i, L in enumerate(fd_layers, 1):
        try:
            t = float(L.thickness)
            gamma = float(L.gamma)
            phi = float(getattr(L, "phi", 0.0))
            # Si la capa es puramente cohesiva (qu>0, phi=0), usar qu/2 = Su.
            qu = float(getattr(L, "qu", 0.0))
            c = float(getattr(L, "c", 0.0))
            if phi == 0 and qu > 0 and c == 0:
                c = qu / 2.0
        except AttributeError as exc:
            raise GeotechError(
                f"Capa #{i} no expone los atributos requeridos: {exc}"
            ) from exc
        if t <= 0:
            raise GeotechError(f"Capa #{i} con thickness <= 0 ({t}).")
        acc += t
        layers.append(SoilLayer(
            name=f"{name_prefix}_{i}",
            cohesion=c, friction_angle=phi, unit_weight=gamma,
        ))
        depths.append(acc)
    return layers, depths


def slope_from_foundation_spec(
    fd_layers: Sequence[Any],
    profile: list[tuple[float, float]],
    *,
    water_table_y: float | None = None,
    surcharge: float = 0.0,
) -> Slope:
    """Construye un ``Slope`` listo para analisis LEM a partir de capas FD.

    Parameters
    ----------
    fd_layers : Sequence
        Capas estratigraficas extraidas del estudio para la fundacion.
    profile : list[tuple]
        Geometria 2D del talud (toe -> ... -> crest).
    water_table_y, surcharge : opcional
        Datos hidraulicos / sobrecarga adicionales.
    """
    layers, depths = layers_from_foundation_spec(fd_layers)
    return Slope(
        profile=profile,
        layers=layers,
        layer_depths=depths,
        water_table_y=water_table_y,
        surcharge=surcharge,
    )


def try_import_fd_bearing() -> Any:
    """Importa ``foundationcalc.services.bearing`` de forma diferida.

    Lanza ``DependencyError`` si el paquete no esta disponible.  Util cuando
    el caller quiere reutilizar el tipo ``SoilLayerSpec`` real (no solo via
    duck typing) y necesita una senal clara.
    """
    try:
        from foundationcalc.services import bearing as _bearing  # type: ignore
        return _bearing
    except ImportError as exc:
        raise DependencyError(
            "foundationcalc no esta instalado. "
            "Instala con: pip install foundationcalc"
        ) from exc


__all__ = [
    "layers_from_foundation_spec",
    "slope_from_foundation_spec",
    "try_import_fd_bearing",
]
