# -*- coding: utf-8 -*-
"""Golden tests — cuña pseudo-estática (mononobe.py) anclada a W00678/W00704.

OJO nomenclatura: ``slipcircle/mononobe.py`` NO implementa el K_AE de
Mononobe-Okabe (empuje sísmico en muros) sino una CUÑA pseudo-estática de
talud indefinido. Los goldens de K_AE de las fichas (K=0.446/0.471,
Ed=154.7/137.5 kN/m) pertenecen a un futuro módulo seismic_foundation y NO
se anclan aquí. Lo que SÍ comparten ambos modelos —y sí está impreso en las
fichas— es el ÁNGULO de rotación sísmica de la gravedad aparente:

    ψ (=θ en las fichas) = arctan( kh / (1 ± kv) )

Fichas máster (fuente real):
  * W00678 — "Seismic design of foundations" (UPC-MEEC, EN 1998-5), pp. de
    empuje sísmico: "Para +kv: tanθ = kh/(1+kv) = 0.2/1.1 → θ=10.3°;
    Para -kv: tanθ = kh/(1-kv) = 0.2/0.9 → θ=12.5°". WorkID: W00678.
  * W00704 — "Seismic design ... Retaining Wall" (Prof. A. Drougkas, UPC):
    mismos dos casos impresos (10.3° / 12.5°). WorkID: W00704.
  Verificación independiente a mano: atan(0.18182)=10.305°,
  atan(0.22222)=12.529° -> los valores impresos 10.3/12.5 SON correctos
  (sin errata).

Mapeo de convenios: mononobe.py define ψ = atan(kh/(1-kv)); el caso "+kv"
de la ficha (denominador 1.1) se obtiene pasando kv=-0.1 y el caso "-kv"
(denominador 0.9) pasando kv=+0.1.

Magnitud del FS: forma cerrada del talud indefinido (literatura clásica,
p.ej. Taylor/Duncan-Wright) y proyección por COMPONENTES del polígono de
fuerzas — dos vías algebraicas independientes de la rotación (ψ) que usa
la librería.
"""
from __future__ import annotations

import math
from math import atan, cos, degrees, radians, sin, tan

import pytest

from slipcircle.mononobe import MononobeInput, evaluate_mononobe

# Caso base: talud indefinido α=30°, H=5 m, γ=18 kN/m³, c=10 kPa, φ=30°
_BASE = dict(slope_height_m=5.0, slope_angle_deg=30.0, gamma_kN_m3=18.0,
             c_kPa=10.0, phi_deg=30.0)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestSeismicAngleW00678:
    """ψ = arctan(kh/(1±kv)) — valores IMPRESOS en W00678/W00704."""

    def test_psi_minus_kv_case_12_5_deg(self):
        """Ficha (caso -kv): kh=0.2, denominador 0.9 -> θ=12.5° impreso.

        A mano: atan(0.2/0.9)=atan(0.22222)=12.529°."""
        res = evaluate_mononobe(MononobeInput(**_BASE, kh=0.2, kv=0.1))
        assert res.psi_deg == pytest.approx(12.529, abs=5e-3)
        assert res.psi_deg == pytest.approx(12.5, abs=0.05)  # valor de ficha

    def test_psi_plus_kv_case_10_3_deg(self):
        """Ficha (caso +kv): kh=0.2, denominador 1.1 -> θ=10.3° impreso.

        A mano: atan(0.2/1.1)=atan(0.18182)=10.305°."""
        res = evaluate_mononobe(MononobeInput(**_BASE, kh=0.2, kv=-0.1))
        assert res.psi_deg == pytest.approx(10.305, abs=5e-3)
        assert res.psi_deg == pytest.approx(10.3, abs=0.05)  # valor de ficha

    def test_psi_zero_without_seismic(self):
        res = evaluate_mononobe(MononobeInput(**_BASE, kh=0.0, kv=0.0))
        assert res.psi_deg == pytest.approx(0.0, abs=1e-12)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestInfiniteSlopeClosedForm:
    """FS estático — forma cerrada clásica del talud indefinido seco:

        FS = c/(γ·H·sinα·cosα) + tanφ/tanα

    A mano (α=30°, H=5, γ=18, c=10, φ=30):
        c-término: 10/(18·5·0.5·0.86603) = 10/38.971 = 0.25660
        φ-término: tan30/tan30 = 1.0
        FS = 1.2566
    """

    def test_static_fs_closed_form(self):
        res = evaluate_mononobe(MononobeInput(**_BASE, kh=0.0, kv=0.0))
        a = radians(30.0)
        fs_closed = (10.0 / (18.0 * 5.0 * sin(a) * cos(a))
                     + tan(radians(30.0)) / tan(a))
        assert fs_closed == pytest.approx(1.2566, abs=2e-4)  # a mano
        assert res.FS_static == pytest.approx(fs_closed, rel=1e-9)
        # sin sismo ambos FS coinciden
        assert res.FS_seismic == pytest.approx(res.FS_static, rel=1e-9)

    def test_pure_friction_slope_fs_tan_ratio(self):
        """c=0: FS = tanφ/tanα exacto (cerrada más simple).
        φ=35°, α=25°: tan35/tan25 = 0.70021/0.46631 = 1.5016."""
        res = evaluate_mononobe(MononobeInput(
            slope_height_m=8.0, slope_angle_deg=25.0, gamma_kN_m3=20.0,
            c_kPa=0.0, phi_deg=35.0, kh=0.0, kv=0.0))
        assert res.FS_static == pytest.approx(1.5016, abs=2e-4)


@pytest.mark.unit
@pytest.mark.safety_critical
class TestPseudoStaticComponentOracle:
    """FS sísmico — oráculo por PROYECCIÓN DE COMPONENTES (independiente de
    la rotación ψ que implementa la librería):

        T_par  = W·[(1-kv)·sinα + kh·cosα]      (motor, paralelo al talud)
        N_perp = W·[(1-kv)·cosα - kh·sinα]      (normal)
        FS     = (c·L + N_perp·tanφ) / T_par ,  L = 1/cosα, W = γ·H

    A mano (caso base, kh=0.2, kv=0.1): W=90, L=1.15470
        T_par  = 90·(0.9·0.5 + 0.2·0.86603) = 56.0885
        N_perp = 90·(0.9·0.86603 - 0.2·0.5) = 61.1481
        FS = (11.5470 + 61.1481·0.57735)/56.0885 = 46.851/56.089 = 0.83530
    """

    def test_seismic_fs_component_projection(self):
        res = evaluate_mononobe(MononobeInput(**_BASE, kh=0.2, kv=0.1))
        a = radians(30.0)
        W = 18.0 * 5.0
        L = 1.0 / cos(a)
        T_par = W * ((1 - 0.1) * sin(a) + 0.2 * cos(a))
        N_perp = W * ((1 - 0.1) * cos(a) - 0.2 * sin(a))
        fs_hand = (10.0 * L + N_perp * tan(radians(30.0))) / T_par
        assert fs_hand == pytest.approx(0.83530, abs=1e-4)  # a mano
        assert res.FS_seismic == pytest.approx(fs_hand, rel=1e-9)
        assert res.is_safe is False  # 0.835 < FS_min=1.1

    def test_component_and_rotation_routes_identical_sweep(self):
        """Identidad algebraica R·sin(α+ψ) = W[(1-kv)sinα + kh·cosα] en un
        barrido de kh/kv: atrapa errores de signo/convención de rotación."""
        for kh in (0.05, 0.15, 0.25):
            for kv in (-0.1, 0.0, 0.1):
                res = evaluate_mononobe(MononobeInput(**_BASE, kh=kh, kv=kv))
                a = radians(30.0)
                W = 18.0 * 5.0
                L = 1.0 / cos(a)
                T_par = W * ((1 - kv) * sin(a) + kh * cos(a))
                N_perp = W * ((1 - kv) * cos(a) - kh * sin(a))
                fs_hand = (10.0 * L + N_perp * tan(radians(30.0))) / T_par
                assert res.FS_seismic == pytest.approx(fs_hand, rel=1e-9), (
                    f"kh={kh} kv={kv}: rotación ψ != proyección componentes"
                )
