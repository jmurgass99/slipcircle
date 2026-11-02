# Changelog

## 1.3.2 — 2026-11-02

- `MorgensternPriceAnalysis`: con f(x) variable el método usaba la resultante Q_i de Spencer con una inclinación θ_i
  distinta en cada dovela, lo que solo impone ΣQ·cos θ = 0: no se cumplía el equilibrio vertical y el FS dependía del
  punto respecto al que se toman momentos. Ahora el equilibrio de fuerzas es dovela a dovela en horizontal y en
  vertical, con X = λ·f(x)·E en cada borde (formulación general de Fredlund y Krahn, 1977): se avanza desde E = 0 en
  el primer borde y el residuo de fuerzas es la E que sobra en el último. f(x) se evalúa en los bordes de las dovelas.
  El momento se toma respecto al centro del círculo o, en una polilínea, respecto al centro del círculo que pasa por
  sus extremos y por su punto medio (con el centroide de la masa el residuo de momentos no es monótono en FS y la
  raíz se pierde). Talud de Duncan y Wright con medio seno: círculo 1,3674 → 1,3587 y polilínea 1,3538 → 1,3585
  (diferencia 1,0 % → 0,015 %); f trapezoidal 1,3638 → 1,3585; f constante sin cambios (≡ Spencer, 1,3587). El test
  `xfail` pasa a ser un test normal y hay otro que comprueba el momento nulo respecto a cuatro puntos en la solución.
  Cambian los casos de uso con M-P: polilínea estratificada 1,5842 → 1,5622 (Spencer 1,5596); talud de Duncan y Wright
  con k_h = 0,10: 1,084 → 1,097 (Spencer 1,098); círculo crítico de Duncan y Wright: 1,3136 → 1,3202 (Spencer 1,3205).

## 1.3.1 — 2026-11-02

- `slices.slices_from_circle` y `Slope.discretize_circle` recortaban el círculo al tramo de perfil definido: un círculo profundo que sale por delante del pie perdía la cuña pasiva (FS = 0,42 frente a 1,11 a mano en un talud φ = 0) y uno que sale más allá del último punto de coronación perdía la parte de detrás. Ahora el terreno se prolonga en horizontal fuera del perfil (`Slope.circle_ground_extent`) y la sobrecarga solo actúa dentro del perfil. Cambian los FS de los casos de uso cuyo círculo pasa por debajo del pie: Duncan y Wright pasa de 1,250 a 1,323 con Bishop (la referencia publicada está en 1,30–1,40).
- `Slope.layer_reference = "elevation"` (opcional): recupera estratos horizontales bajo y = 0 para un talud que corta estratos horizontales, como un terraplén sobre arcilla blanda; el valor por defecto `"surface"` mide desde la superficie local.
- Casos de uso: el freático ya no sube por encima del pie (el paquete no aplica el peso ni los empujes del agua exterior); el ejemplo del terraplén usa `layer_reference="elevation"` y compara con N_c·c_u/(γ·H).
- Test `xfail` estricto que documenta la diferencia círculo/polilínea de M-P con f de medio seno (≈ 1 %): la formulación Q impone el equilibrio horizontal y de momentos pero no ΣQ·sen θ = 0.

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
