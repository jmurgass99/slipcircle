# Changelog

## 1.3.0 — 2026-10-31

Correcciones encontradas al ejecutar los casos de uso del README.

- `spencer.SpencerAnalysis` y `morgenstern_price.MorgensternPriceAnalysis`: el numerador de la resultante entre dovelas Q_i era el de Bishop (c·b y W − u·b) y no el de Spencer (1967) Ec. 7 (c·l y W·cos α − u·l con l = b·sec α); en φ = 0 no reproducían la forma cerrada c·R·L_arco/ΣW·(x_m − c_x). Además la iteración alternaba raíces en FS y en θ (o λ) con un bracket de FS que incluía los polos de Q_i, de modo que el resultado dependía del número de dovelas (θ se quedaba en 0 sin equilibrio de momentos o se disparaba a FS = 8,7 con θ = −46°) y `converged` era siempre False. Ahora se usa el esquema clásico FS_f(θ) = FS_m(θ) (`_parallel_solver`) con el bracket de FS por encima de los polos. En el talud de Duncan y Wright ambos convergen a 1,245–1,250 (Bishop 1,250) con 30, 40 y 60 dovelas; en φ = 0 dan 1,2570, la forma cerrada; M-P con f constante coincide con Spencer. El resultado expone `FS_force_eq` y `FS_moment_eq`. La rama de polilínea sumaba el momento de las fuerzas entre dovelas (interiores) y mezclaba convenios de signo entre W y N/S: una polilínea muestreada del círculo de Duncan y Wright no reproducía su FS; ahora el momento global de la masa (peso, sísmica y reacciones de la base) es común a las dos ramas y da 1,243 frente a 1,245. Con acción sísmica, W_h entraba con signo opuesto en N/S respecto a Q y el método no convergía con k_h ≠ 0; corregido (Spencer con k_h = 0,10 converge a pocos % del Bishop pseudo-estático).
- `api.app`: el endpoint `POST /api/slope/newmark` existía en `endpoints_newmark.py` pero nunca se registraba en la aplicación (404); registrado.
- Hallazgo documentado, sin cambio de código: `seismic.PseudoStaticAnalysis` (Bishop pseudo-estático) aplica W_h en la base de la dovela (R·k_h·W·cos α) y usa W sin el factor (1 − k_v) en el momento motriz, mientras que la resistencia sí lleva W·(1 − k_v). Con k_h = 0,10 y k_v = 0,05 da 0,956; Spencer (W_h en el centroide) 1,037; Spencer con W_h en la base 1,004 y la fórmula documentada de Bishop 1,005. Es el método más conservador de los tres y se mantiene por compatibilidad; la diferencia queda recogida en el README.
- `janbu.JanbuAnalysis` acepta superficies polilínea (antes lanzaba NotImplementedError pese a anunciarlas); el refuerzo solo se contabiliza con círculos.
- 25 casos de uso ejecutados en `examples/` con su salida en `examples/output/`; README reescrito; metadatos del paquete alineados.

## 1.2.0

Goldens de la campaña del máster: forma cerrada no drenada, talud infinito, pulso de Newmark, cuña pseudo-estática (W00678/W00704), geotextil (W00222), Bray-Travasarou con intercepto por variante y bloque de Newmark sin «planeo».

## 1.0.0

Publicación inicial como repositorio independiente (origen: SlopeStabilityJm del monorepo Structural JM).
