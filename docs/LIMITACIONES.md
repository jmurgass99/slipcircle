# Limitaciones conocidas de slipcircle

Lo que el paquete no hace o hace de forma simplificada.

- La búsqueda del círculo crítico es un barrido en malla de centros y radios (Bishop por defecto; `method=` admite Spencer, Morgenstern-Price y Janbu) con un refinado local Nelder-Mead opcional (`optimize_circle`); las superficies no circulares no entran en la búsqueda.
- `Slope.layer_at_depth` mide por defecto la profundidad desde la superficie local (`layer_reference="surface"`): las interfaces son paralelas al perfil. Para estratos horizontales cortados por el talud (terraplén sobre blando, desmonte en estratos horizontales) usa `layer_reference="elevation"` con las profundidades bajo y = 0. `slope_from_foundation_spec` crea el talud con el valor por defecto.
- Las presiones de poros son datos (freático, piezométrica, r_u); no hay red de flujo. No se modela agua exterior: un freático por encima del terreno delante del pie da presiones intersticiales sin el peso ni los empujes de la lámina de agua y un FS artificialmente bajo.
- `PseudoStaticAnalysis` (Bishop pseudo-estático) aplica W_h en la base de la dovela y usa W sin (1 − k_v) en el momento motriz: es el más conservador. Spencer y M-P con `kh` aplican W_h en el centroide. Las cifras en el talud de Duncan y Wright con k_h = 0,10 y k_v = 0,05 están en la última línea de su caso de uso.
- La estimación de desplazamiento de `PseudoStaticAnalysis.newmark_displacement` (Ambraseys-Menu) usa k_y = k_h del análisis, no la aceleración de fluencia del talud; es de cribado.
- `mononobe` es una cuña de talud indefinido con la gravedad aparente girada ψ, no el K_AE de un muro.
- El refuerzo entra como momento resistente constante (tensión de diseño limitada por anclaje o por el elemento), sin interacción con la deformación; con polilíneas solo lo aplica Janbu generalizado.
- `drains_casagrande` y `fs_by_code` son fórmulas y tablas de anteproyecto.
- `bray_travasarou` se declara `experimental`: ecuación publicada reproducida, sin calibración propia.
