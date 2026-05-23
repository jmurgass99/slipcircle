# SlopeStabilityJm — Estabilidad de Taludes Multi-Método

[![CI](https://github.com/jmurgass99/slopestabilityjm/actions/workflows/ci.yml/badge.svg)](https://github.com/jmurgass99/slopestabilityjm/actions/workflows/ci.yml)
[![Coverage](https://img.shields.io/badge/coverage-84.31%25-brightgreen?style=flat)](https://github.com/jmurgass99/slopestabilityjm)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-Proprietary-red)]()

## Objetivo
Librería Python para el análisis de estabilidad de taludes mediante métodos de equilibrio límite (LEMs). Implementa los métodos de **Bishop Simplificado**, **Janbu Generalizado**, **Morgenstern-Price** y **Spencer**, con soporte para análisis pseudo-estático sísmico y refuerzo con geosintéticos.

Este proyecto llena un **vacío geotécnico** del ecosistema: `FoundationDesignJm` cubre muros de contención (Rankine, Coulomb, Mononobe-Okabe) pero NO analiza taludes naturales ni cortes.

---

## ¿Qué va a hacer?

### 1. Definición del Talud
- Geometría por puntos (perfil del terreno, capas de suelo)
- Múltiples estratos con parámetros c', φ', γ distintos
- Nivel freático: posición fija, línea piezométrica o presión de poros ru
- Sobrecargas en corona y berma

### 2. Métodos de Análisis (Equilibrio Límite)
- **Bishop Simplificado:** Superficies circulares, satisface equilibrio de momentos. Iterativo en FS.
- **Janbu Generalizado:** Superficies no circulares, satisface equilibrio de fuerzas.
- **Spencer:** Satisface equilibrio de momentos Y fuerzas. Asume inclinación constante de fuerzas entre dovelas.
- **Morgenstern-Price:** Generalización de Spencer con función f(x) arbitraria para la inclinación de fuerzas entre dovelas. El más riguroso.

### 3. Búsqueda de Superficie Crítica
- Barrido de centros y radios (superficies circulares)
- Algoritmo de optimización para la superficie con FS mínimo
- Grid search + refinamiento con gradiente

### 4. Análisis Sísmico (Pseudo-Estático)
- Coeficientes sísmicos horizontales y verticales: kh, kv
- Recomendaciones EC8-5: kh = 0.5·ag/g (o según desplazamiento admisible)
- Método de Newmark (desplazamiento permanente por integración doble)

### 5. Refuerzo
- Geomallas / geotextiles: fuerza de refuerzo Tr por capa
- Anclajes (soil nails): fuerza axial + inclinación
- Verificación de mecanismos de falla: pull-out, rotura, deslizamiento sobre refuerzo

### 6. Visualización
- Renderizado del perfil del talud con capas de suelo coloreadas
- Superficie de falla crítica superpuesta
- Distribución de fuerzas por dovela

---

## Estructura de Módulos Planificada

```
slopestabilityjm/
├── __init__.py
├── geometry.py             # Perfil del talud, capas, nivel freático
├── soil_profile.py         # Estratos: c', φ', γ, γsat, ru
├── slices.py               # Discretización en dovelas (slice method)
├── bishop.py               # Bishop simplificado (circular, momentos)
├── janbu.py                # Janbu generalizado (no circular, fuerzas)
├── spencer.py              # Spencer (circular, momentos + fuerzas)
├── morgenstern_price.py    # Morgenstern-Price (general, f(x) arbitraria)
├── search.py               # Búsqueda de superficie crítica (grid + optimización)
├── seismic.py              # Pseudo-estático (kh, kv) + Newmark displacement
├── reinforcement.py        # Geomallas, anclajes, soil nails
├── water.py                # Presión de poros: nivel freático, ru, piezométrica
├── visualization.py        # Renderizado matplotlib del talud + falla + dovelas
└── tests/
    ├── test_bishop.py
    ├── test_janbu.py
    ├── test_spencer.py
    ├── test_seismic.py
    └── test_reinforced.py
```

---

## Fuentes de Información

### Teoría y Referencias Académicas
| Fuente | Contenido |
|--------|-----------|
| **Duncan, Wright & Brandon (2014)** — *Soil Strength and Slope Stability* | Teoría completa de LEMs, Bishop, Spencer, Morgenstern-Price |
| **Abramson et al. (2002)** — *Slope Stability and Stabilization Methods* | Refuerzo, anclajes, geosintéticos |
| **EC7 EN 1997-1** | Factores parciales para estabilidad de taludes (DA-1, DA-3) |
| **EC8-5 EN 1998-5** | Análisis pseudo-estático de taludes, kh y desplazamiento Newmark |
| **Curso UPC — Foundation Structures** | Empujes de tierras, Mononobe-Okabe (ya en `informacion para modelar`) |

### Repositorios de Referencia (ya descargados en `Nuevos/`)
- `PySlope-main` — Implementación en Python del método de Bishop. **Código fuente principal de referencia.** Usa superficies circulares con grid search.
- `FoundationDesign-main` (en `Usados/`) — Implementación de Rankine/Coulomb para muros, referencia complementaria.

### Integración con FoundationDesignJm
- `foundationdesignjm/core/` — Modelo `SoilLayer` y `SoilProfile` reutilizable
- `foundationdesignjm/analysis/earth_pressure.py` — Parámetros de suelo ya definidos
- `foundationdesignjm/analysis/seismic.py` — Mononobe-Okabe ya implementado

### Integración con DynamicsJmPro (Web)
- Futuro: módulo de visualización web del talud con Three.js (2D/3D)
- Análisis paramétrico: variar kh y ver cómo cambia FS en tiempo real

---

## Criterio de Éxito
- [ ] Calcular FS de un talud homogéneo con Bishop y obtener resultado dentro de ±2% del ejemplo de Duncan & Wright
- [ ] Manejar taludes multicapa con nivel freático
- [ ] Análisis pseudo-estático con kh según EC8-5
- [ ] Búsqueda automática de la superficie crítica
- [ ] Visualización del talud con la superficie de falla y FS
- [ ] Suite de tests con ≥25 casos validados contra literatura

## Ecosystem

This project is part of the **jmurgass99 structural engineering ecosystem**:

- [fastsap](https://github.com/jmurgass99/fastsap) - orquestador SAP2000/ETABS
- [dynamicsjmpro](https://github.com/jmurgass99/dynamicsjmpro) - SaaS frontend + backend
- [dynamicsjm](https://github.com/jmurgass99/dynamicsjm) - motor dinamico
- [foundationdesignjm](https://github.com/jmurgass99/foundationdesignjm) - cimentaciones
- [windloadjm](https://github.com/jmurgass99/windloadjm) - cargas de viento
- [seismicrecordsjm](https://github.com/jmurgass99/seismicrecordsjm) - records sismicos
- [slopestabilityjm](https://github.com/jmurgass99/slopestabilityjm) - estabilidad de taludes
- [steelconnectionjm](https://github.com/jmurgass99/steelconnectionjm) - conexiones de acero
- [steeloptimizerjm](https://github.com/jmurgass99/steeloptimizerjm) - optimizacion GA
- [coldformedjm](https://github.com/jmurgass99/coldformedjm) - acero conformado en frio
- [fatiguejm](https://github.com/jmurgass99/fatiguejm) - analisis de fatiga
- [firesteeljm](https://github.com/jmurgass99/firesteeljm) - analisis termico al fuego
- [dxfstructuraljm](https://github.com/jmurgass99/dxfstructuraljm) - generacion DXF
