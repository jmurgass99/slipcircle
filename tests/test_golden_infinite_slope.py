# -*- coding: utf-8 -*-
"""Golden de forma cerrada: TALUD INFINITO SECO.

Motor bajo test
---------------
`slices.slices_from_polyline` (discretizacion: pesos, alturas, inclinacion de
base) + `fellenius.fellenius_fs` (reductor Ordinario de Dovelas). Es la unica
via del paquete que admite una superficie de rotura PLANA: `BishopAnalysis` y
`JanbuAnalysis` solo aceptan `CircularSurface`, con la que el caso de talud
infinito no es representable. Se deja constancia de esa limitacion aqui.

Procedencia de los valores esperados
------------------------------------
GOLDEN. Formas cerradas de talud infinito, derivadas abajo a mano y aseveradas
literalmente. NINGUNA cifra procede de ejecutar la implementacion.

Derivacion (talud infinito seco)
--------------------------------
Superficie de rotura plana paralela a la superficie del terreno, inclinada un
angulo beta, a profundidad VERTICAL z. Todas las dovelas son identicas, y por
periodicidad las fuerzas entre dovelas se cancelan: el equilibrio limite es
exacto (no aproximado) para este caso.

Dovela de ancho b:
    W = gamma * z * b            (peso, altura vertical z)
    l = b / cos(beta)            (longitud de la base)
    N = W * cos(beta)            (normal a la base, seco -> u = 0)
    T_motriz = W * sin(beta)

    FS = (c*l + N*tan(phi)) / (W*sin(beta))
       = (c*b/cos(beta) + gamma*z*b*cos(beta)*tan(phi)) / (gamma*z*b*sin(beta))

Multiplicando numerador y denominador por cos(beta) se recupera la forma
clasica de manual:

    FS = (c + gamma*z*cos^2(beta)*tan(phi)) / (gamma*z*cos(beta)*sin(beta))

y con c = 0 el peso, el ancho y la profundidad se cancelan por completo:

    FS = tan(phi) / tan(beta)        <- caso limite exigido

Numeros elegidos (beta = 15 deg, phi = 30 deg, gamma = 20 kN/m3, z = 3 m)
------------------------------------------------------------------------
tan(30 deg) = 1/sqrt(3);  tan(15 deg) = 2 - sqrt(3)

  (a) c' = 0:
      FS = (1/sqrt(3)) / (2 - sqrt(3)) = 1 / (2*sqrt(3) - 3) = 1 + 2/sqrt(3)
         = 2.1547005383792515

  (b) c' = 10 kPa:
      cos^2(15) = (1 + cos 30)/2 = (2 + sqrt(3))/4
      cos(15)*sin(15) = sin(30)/2 = 1/4
      gamma*z = 60 kPa
      numerador   = 10 + 60*(2+sqrt(3))/4 * (1/sqrt(3)) = 10 + 15*(2/sqrt(3) + 1)
                  = 25 + 10*sqrt(3) = 42.320508075688775
      denominador = 60/4 = 15
      FS = (25 + 10*sqrt(3))/15 = (5 + 2*sqrt(3))/3 = 2.821367205045918

Tolerancia
----------
rel=1e-12: el resultado es EXACTO (aritmetica de punto flotante aparte). No hay
error de discretizacion porque todas las dovelas tienen la misma alpha = beta y
la misma altura, luego el numero de dovelas no interviene; se comprueba
explicitamente en `test_fs_independent_of_discretisation_and_depth`.
"""
from __future__ import annotations

import math

import pytest

from slopestabilityjm.fellenius import Slice, fellenius_fs
from slopestabilityjm.geometry import Slope
from slopestabilityjm.slices import slices_from_polyline
from slopestabilityjm.soil_profile import SoilLayer

_BETA_DEG = 15.0
_PHI_DEG = 30.0
_GAMMA = 20.0          # kN/m3
_Z = 3.0               # profundidad VERTICAL del plano de rotura [m]
_TAN_BETA = math.tan(math.radians(_BETA_DEG))

# Tramo central del talud usado para el analisis (el resultado no depende de el)
_X0, _X1 = 20.0, 80.0


def _infinite_slope_slices(c_kPa: float, n_slices: int = 40,
                           z: float = _Z) -> list[Slice]:
    """Ladera recta de pendiente beta con plano de rotura paralelo a z metros.

    El perfil del terreno es la recta y = x*tan(beta) entre x=0 y x=100; la
    superficie de rotura es la recta paralela y = x*tan(beta) - z. Al ser ambas
    rectas paralelas, cada dovela tiene altura z y base inclinada beta: es la
    definicion de talud infinito.
    """
    layer = SoilLayer("Talud infinito", cohesion=c_kPa,
                      friction_angle=_PHI_DEG, unit_weight=_GAMMA)
    slope = Slope(
        profile=[(0.0, 0.0), (100.0, 100.0 * _TAN_BETA)],
        layers=[layer],
        layer_depths=[200.0],   # un unico estrato, mas profundo que el corte
        water_table_y=None,     # SECO -> u = 0
    )
    polyline = [(_X0, _X0 * _TAN_BETA - z), (_X1, _X1 * _TAN_BETA - z)]
    dicts = slices_from_polyline(slope, polyline, n_slices=n_slices)
    assert dicts, "la discretizacion no produjo dovelas"
    return [
        Slice(
            W_kN=d["weight_kN"],
            alpha_deg=math.degrees(d["alpha_rad"]),
            c_kPa=d["cohesion"],
            phi_deg=math.degrees(d["phi_rad"]),
            base_length_m=d["width"] / math.cos(d["alpha_rad"]),
            u_kPa=d["u_kPa"],
        )
        for d in dicts
    ]


class TestInfiniteSlopeDry:
    """FS = tan(phi)/tan(beta) para talud infinito seco sin cohesion."""

    def test_geometry_is_actually_an_infinite_slope(self):
        """Precondicion: toda dovela tiene alpha = beta = 15 deg, h = z = 3 m
        y u = 0. Si esto falla, el golden de FS no estaria midiendo lo que dice.
        """
        slices = _infinite_slope_slices(c_kPa=0.0)
        for s in slices:
            assert s.alpha_deg == pytest.approx(_BETA_DEG, rel=1e-12)
            assert s.u_kPa == 0.0
        # W = gamma*z*b, con b = (80-20)/40 = 1.5 m -> 20*3*1.5 = 90 kN/m
        assert slices[0].W_kN == pytest.approx(_GAMMA * _Z * 1.5, rel=1e-12)
        assert slices[0].W_kN == pytest.approx(90.0, rel=1e-12)

    def test_fs_equals_tanphi_over_tanbeta(self):
        """FS = tan(30)/tan(15) = 1 + 2/sqrt(3) = 2.1547005383792515."""
        fs_closed = math.tan(math.radians(_PHI_DEG)) / _TAN_BETA
        # la forma cerrada coincide con el numero derivado a mano
        assert fs_closed == pytest.approx(1.0 + 2.0 / math.sqrt(3.0), rel=1e-12)
        assert fs_closed == pytest.approx(2.1547005383792515, rel=1e-12)

        fs = fellenius_fs(_infinite_slope_slices(c_kPa=0.0))
        assert fs == pytest.approx(2.1547005383792515, rel=1e-12)

    def test_fs_independent_of_discretisation_and_depth(self):
        """En c'=0 el FS no depende ni del numero de dovelas ni de z.

        Justifica la tolerancia rel=1e-12: no hay error de discretizacion que
        tolerar, porque z y b se cancelan analiticamente en el cociente.
        """
        ref = 2.1547005383792515
        for n in (3, 7, 40, 250):
            assert fellenius_fs(_infinite_slope_slices(0.0, n_slices=n)) == \
                pytest.approx(ref, rel=1e-12)
        for z in (0.5, 1.5, 3.0, 8.0):
            assert fellenius_fs(_infinite_slope_slices(0.0, z=z)) == \
                pytest.approx(ref, rel=1e-12)

    def test_fs_with_cohesion_closed_form(self):
        """c'=10 kPa: FS = (5 + 2*sqrt(3))/3 = 2.821367205045918.

        Mismo caso con cohesion, para que el golden no dependa solo del termino
        friccional (c'=0 cancela peso, ancho y profundidad; aqui NO se cancelan
        y el termino cohesivo ejercita los pesos y anchos que produce
        `slices_from_polyline`).

        ALCANCE HONESTO: la longitud de base ``l = b/cos(alpha)`` la calcula
        ESTE test, no el paquete — `slices_from_polyline` no emite `l` y
        `fellenius.Slice` la exige. Los que sí proceden del paquete son
        ``b = width``, ``alpha_rad``, ``W = weight_kN`` y ``u``. Es decir, este
        golden ancla la discretizacion y el reductor, no la formula de l.
        """
        beta = math.radians(_BETA_DEG)
        fs_closed = (
            (10.0 + _GAMMA * _Z * math.cos(beta) ** 2
             * math.tan(math.radians(_PHI_DEG)))
            / (_GAMMA * _Z * math.cos(beta) * math.sin(beta))
        )
        assert fs_closed == pytest.approx((5.0 + 2.0 * math.sqrt(3.0)) / 3.0,
                                          rel=1e-12)
        assert fs_closed == pytest.approx(2.821367205045918, rel=1e-12)

        fs = fellenius_fs(_infinite_slope_slices(c_kPa=10.0))
        assert fs == pytest.approx(2.821367205045918, rel=1e-12)

    def test_limit_state_phi_equals_beta(self):
        """Caso limite estricto: phi = beta -> FS = 1 exacto (talud en reposo
        justo en el angulo de rozamiento interno). Se construye con beta=30 y
        phi=30 reutilizando la misma derivacion.
        """
        layer = SoilLayer("Arena", cohesion=0.0, friction_angle=30.0,
                          unit_weight=19.0)
        t30 = math.tan(math.radians(30.0))
        slope = Slope(profile=[(0.0, 0.0), (100.0, 100.0 * t30)],
                      layers=[layer], layer_depths=[200.0])
        poly = [(20.0, 20.0 * t30 - 2.0), (80.0, 80.0 * t30 - 2.0)]
        ds = slices_from_polyline(slope, poly, n_slices=30)
        slices = [
            Slice(W_kN=d["weight_kN"], alpha_deg=math.degrees(d["alpha_rad"]),
                  c_kPa=d["cohesion"], phi_deg=math.degrees(d["phi_rad"]),
                  base_length_m=d["width"] / math.cos(d["alpha_rad"]),
                  u_kPa=d["u_kPa"])
            for d in ds
        ]
        assert fellenius_fs(slices) == pytest.approx(1.0, rel=1e-12)
