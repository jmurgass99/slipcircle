# Ejecuta la suite completa con cobertura y aplica el gate (fail_under en pyproject).
# Uso:  .\scripts\run_coverage.ps1
$ErrorActionPreference = "Stop"
Push-Location $PSScriptRoot\..
try {
    python -m pytest tests/ --cov=slopestabilityjm --cov-report=term-missing --cov-report=html
    Write-Host ""
    Write-Host "Reporte HTML generado en: htmlcov/index.html"
} finally {
    Pop-Location
}
