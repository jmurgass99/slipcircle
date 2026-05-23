# -*- coding: utf-8 -*-
"""Metodo Ordinario de Dovelas (Fellenius)."""
import math

import pytest

from slopestabilityjm.fellenius import Slice, fellenius_fs


def _slices():
    # Talud simple, 4 dovelas (datos ilustrativos, suelo c-phi).
    return [
        Slice(W_kN=120.0, alpha_deg=-5.0, c_kPa=10.0, phi_deg=28.0, base_length_m=2.2),
        Slice(W_kN=180.0, alpha_deg=10.0, c_kPa=10.0, phi_deg=28.0, base_length_m=2.0),
        Slice(W_kN=160.0, alpha_deg=25.0, c_kPa=10.0, phi_deg=28.0, base_length_m=2.1),
        Slice(W_kN=90.0,  alpha_deg=40.0, c_kPa=10.0, phi_deg=28.0, base_length_m=2.6),
    ]


class TestFellenius:
    def test_fs_positivo_y_formula(self):
        sl = _slices()
        fs = fellenius_fs(sl)
        # Recalculo manual de la formula.
        R = sum(s.c_kPa * s.base_length_m
                + max(s.W_kN * math.cos(math.radians(s.alpha_deg))
                      - s.u_kPa * s.base_length_m, 0.0) * math.tan(math.radians(s.phi_deg))
                for s in sl)
        D = sum(s.W_kN * math.sin(math.radians(s.alpha_deg)) for s in sl)
        assert fs == pytest.approx(R / D, rel=1e-9)
        assert fs > 0

    def test_presion_poros_reduce_fs(self):
        seco = fellenius_fs(_slices())
        con_agua = fellenius_fs([Slice(s.W_kN, s.alpha_deg, s.c_kPa, s.phi_deg,
                                       s.base_length_m, u_kPa=20.0) for s in _slices()])
        assert con_agua < seco

    def test_sismo_reduce_fs(self):
        estatico = fellenius_fs(_slices())
        sismico = fellenius_fs(_slices(), seismic_kh=0.15)
        assert sismico < estatico

    def test_sin_dovelas_lanza(self):
        with pytest.raises(ValueError):
            fellenius_fs([])
