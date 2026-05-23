"""Tests de la jerarquia de excepciones tipadas (Mejora B)."""
from __future__ import annotations
import pytest

from slopestabilityjm.exceptions import (
    SlopeError, GeotechError, ConvergenceError, CodeError,
    GeometryError, ConfigurationError, DependencyError,
)


class TestExceptionHierarchy:

    def test_all_inherit_from_slope_error(self):
        for cls in (GeotechError, ConvergenceError, CodeError,
                    GeometryError, ConfigurationError, DependencyError):
            assert issubclass(cls, SlopeError)

    def test_slope_error_is_value_error(self):
        """Backward compat: codigo legacy que captura ValueError sigue funcionando."""
        assert issubclass(SlopeError, ValueError)
        for cls in (GeotechError, ConvergenceError, GeometryError):
            assert issubclass(cls, ValueError)

    def test_convergence_error_carries_metadata(self):
        exc = ConvergenceError("no convergio", iteration=42, residual=1.23e-3)
        assert exc.iteration == 42
        assert exc.residual == pytest.approx(1.23e-3)
        assert "no convergio" in str(exc)

    def test_convergence_error_defaults(self):
        exc = ConvergenceError("x")
        assert exc.iteration == -1
        import math
        assert math.isnan(exc.residual)

    def test_can_catch_via_value_error(self):
        with pytest.raises(ValueError):
            raise GeotechError("phi negativo")

    def test_can_catch_via_specific_type(self):
        with pytest.raises(GeotechError):
            raise GeotechError("phi negativo")


class TestBishopStrictMode:
    """Bishop en modo strict lanza ConvergenceError cuando max_iter se agota."""

    def test_strict_mode_raises_on_no_convergence(self, simple_slope, test_surface):
        from slopestabilityjm.bishop import BishopAnalysis
        # max_iter=1 garantiza no convergencia
        ba = BishopAnalysis(simple_slope, n_slices=30, tol=1e-12,
                            max_iter=1, strict=True)
        with pytest.raises(ConvergenceError) as info:
            ba.factor_of_safety(test_surface)
        assert info.value.iteration == 1
        assert info.value.residual > 0

    def test_non_strict_returns_estimate(self, simple_slope, test_surface):
        from slopestabilityjm.bishop import BishopAnalysis
        ba = BishopAnalysis(simple_slope, n_slices=30, tol=1e-12,
                            max_iter=1, strict=False)
        r = ba.factor_of_safety(test_surface)
        # No raise, devuelve FS estimado
        assert r["FS"] is not None


class TestReExports:
    """Las excepciones se reexportan desde el paquete raiz."""

    def test_root_reexports(self):
        import slopestabilityjm as ss
        assert ss.SlopeError is SlopeError
        assert ss.ConvergenceError is ConvergenceError
        assert ss.GeotechError is GeotechError
        assert ss.GeometryError is GeometryError
        assert ss.CodeError is CodeError
        assert ss.DependencyError is DependencyError
        assert ss.ConfigurationError is ConfigurationError
