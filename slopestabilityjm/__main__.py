"""CLI minimo de SlopeStabilityJm.

Uso:
    slopestabilityjm --version
    slopestabilityjm api        # arranca FastAPI con uvicorn (requiere extra [api])
"""
from __future__ import annotations
import sys
from . import __version__


def main(argv: list[str] | None = None) -> int:
    """Punto de entrada del script ``slopestabilityjm``."""
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help"):
        print("SlopeStabilityJm v{}  --  uso: slopestabilityjm [--version|api]"
              .format(__version__))
        return 0
    if args[0] == "--version":
        print(__version__)
        return 0
    if args[0] == "api":
        try:
            import uvicorn  # type: ignore
            from .api.app import app  # noqa: F401
        except ImportError as exc:
            from .exceptions import DependencyError
            raise DependencyError(
                "FastAPI/uvicorn no instalados. Instala con: "
                "pip install slopestabilityjm[api]"
            ) from exc
        host = "127.0.0.1"
        port = 8000
        uvicorn.run("slopestabilityjm.api.app:app", host=host, port=port)
        return 0
    print(f"Comando desconocido: {args[0]}", file=sys.stderr)
    return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
