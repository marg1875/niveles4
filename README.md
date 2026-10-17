# niveles4 — Caracterización y clasificación de niveles de KMI (MATLAB + Weka)

Migración del pipeline en Python [`marg1875/Niveles3`](https://github.com/marg1875/Niveles3) a una
arquitectura **MATLAB (preprocesamiento + exponentes de Hurst) → ARFF → Weka (clasificación)**,
siguiendo la metodología de:

> Martinez-Peon D., Garcia-Hernandez N. V., Benavides-Bravo F. G., Parra-Vega V. (2024).
> *Characterization and classification of kinesthetic motor imagery levels.* J. Neural Eng. 21, 046024.

Datos: 5 pacientes con ECV isquémico (Px006–Px010), meses 1/3/6, imaginería motora kinestésica
(`*img*`, **KMI**) y ejecución física (**FIS**) al 10 % y 40 % de la FVM, registrados con un
**g.Nautilus de 16 canales a 250 Hz** (`Datos/`).

## Requisitos

- MATLAB (probado en R2023b) con **Signal Processing Toolbox**.
- Opcionales: Parallel Computing Toolbox (paraleliza las etapas 1, 2 y 4), Statistics Toolbox
  (`signrank` exacto; sin él se usa la aproximación normal) y Wavelet Toolbox (una prueba con fGn).
- **Weka 3.8.x** (probado con 3.8.6). Se autodetecta en `C:\Program Files\Weka-*\weka.jar`; si no,
  defina la variable de entorno `WEKA_JAR`. Weka se usa a través de su API Java dentro de MATLAB,
  sin parsear texto, y da los mismos resultados que Weka Explorer o la CLI con `-x 5 -s 1`.

## Uso

```matlab
>> ejecutar_pipeline          % corrida completa (reanudable)
>> qa_alineacion_estimulos    % evidencia de la alineación stims <-> señal
>> runtests('pruebas')        % pruebas unitarias
```

Todo se configura en `config_niveles.m`. En `ejecutar_pipeline.m` hay ajustes rápidos comentados,
por ejemplo `cfg.modo_prueba = true` (solo Px006 mes 1, 2 duraciones y 2 configuraciones; los
resultados van a `resultados_prueba/`). Cada etapa omite lo ya calculado. Para recalcular una etapa,
borre su carpeta en `resultados/`.

| Etapa | Archivo | Salida |
|---|---|---|
| 1 | `etapa1_preprocesar.m` — inventario, filtro, z-score, blanqueo, ICA por sesión | `registros.csv`, `sesiones.csv`, `senales/` |
| 2 | `etapa2_caracteristicas.m` — segmentación multi-duración, HO/HRS(p)/HV | `caracteristicas/caracteristicas_d<N>.csv` |
| 3 | `etapa3_exportar_arff.m` — datasets para Weka | `arff/d<N>/{intra,mes,global}/*.arff` + `*_meta.csv` |
| 4 | `etapa4_clasificar_weka.m` — CV 5-fold en Weka (fases A y B + chequeo de fuga) | `weka/metricas.csv`, `weka/mejor_duracion.txt` |
| 5 | `etapa5_reportes.m` — tablas y figuras | `reportes/*.csv`, `reportes/figuras/*.png` |

Los ARFF abren directamente en **Weka Explorer / Experimenter**. Cada uno contiene los 16 canales
× 9 características (`HO_<canal>`, `HRS<p>_<canal>`, `HV_<canal>`) y la clase `nivel` ∈ {N00, N10, N40}.

## Metodología

1. **Preprocesamiento** (por archivo): quitar DC → notch 60 Hz (Q = 30) → Butterworth pasabanda
   0.5–45 Hz de orden 4 en SOS, fase cero (opción `cfg.filtro.tipo = 'fir'`) → [ASR, opcional] →
   referencia promedio común (CAR) → z-score por canal (Ec. 1).

   **Artefactos.** Los filtros no eliminan los picos: un pico es de banda ancha y cae dentro de
   0.5–45 Hz. Además, el R/S de Hurst es muy sensible a ellos porque mide el *rango* de la ventana.
   - **ASR** (Artifact Subspace Reconstruction; Mullen et al., 2015, IEEE TBME 62:2553; Chang et
     al., 2020, IEEE TBME 67:1114). Implementado en `asr_calibrar.m` y `asr_procesar.m`:
     - Se calibra con las ventanas limpias de cada registro.
     - En ventanas de 0.5 s, las componentes cuya varianza supera k desviaciones de la calibración
       (k = 20; `cfg.asr`) se reconstruyen a partir del resto. **No se descartan ventanas.**
     - Simplificaciones frente a clean_rawdata: no se usa el prefiltro espectral, porque la señal ya
       está en 0.5–45 Hz, y la distribución se resume con mediana y MAD.
     - `comparar_asr.m` compara la señal sin ASR, ASR con k = 20 y ASR con k = 30.
     - **Preprocesamiento final: ASR con k = 30** (`cfg.asr`). Se eligió por ser el valor más
       conservador del rango recomendado: repara 4.3 % del registro, frente a 6.2 % con k = 20.
     - Con k = 30, el LOSO KMI 10 vs 40 es de 59.1 % (p = 0.034 por permutación, 1000
       permutaciones). Con k = 20 da 55.8 % y sin ASR 56.4 % (p = 0.107).
     - El análisis completo con este preprocesamiento está en `corrida_asr.m` → `resultados_asr/`.
   - **Rechazo de ventanas por 6 MAD** (`marcar_artefactos.m`):
     - Criterio: |x − mediana| / (1.4826·MAD) > 6 en algún canal, con margen de ±0.25 s.
     - **Solo como análisis de sensibilidad** (`cfg.artefactos.activar = false`): descarta ~50 % de
       las ventanas y cambia la generalización entre sesiones en apenas ~2 puntos (56.4 % → 58.4 %;
       `comparar_variantes.m`).
   - **Electrodos puenteados** (r > 0.99 tras CAR, concentrados en el mes 6): no se corrigen ni se
     excluyen sesiones. Interpolar, eliminar o excluir no mejoró la generalización
     (`comparar_variantes.m`), así que se reportan como limitación.
2. **Blanqueo + ICA por sesión** (paciente × mes × condición, con los archivos de 10 % y 40 %
   concatenados). Se usa **una W por sesión** para que reposo, 10 % y 40 % pasen por el mismo
   desmezclador.
   - Blanqueo PCA (Ec. 2; opción `p0_articulo` para la Ec. 3).
   - ICA de **gradiente natural con no linealidad flexible** según el signo de la curtosis
     (Ecs. 5–6; α = 4, β = 2), en su forma no holonómica, con W₀ = I (determinista) y paso adaptativo.
   - Se anulan hasta 2 componentes con curtosis atípica (z robusto de log(κ+3) > 2) y se reproyecta
     a electrodos ("sensores limpios").
3. **Segmentación.** Las 80 entradas de `stims` se reparten uniformemente sobre `EEG_all` y forman
   20 bloques de reposo/tarea de ~7.8 s. Tras márgenes de 0.5 s / 0.25 s se fijan **3 anclas por
   bloque**. Para cada duración (512 muestras, como el artículo, y 2, 1.5, 1, 0.5, 0.4, 0.3, 0.2 y
   0.1 s) se toma la ventana centrada en cada ancla: mismos instantes y mismo número de instancias en
   todas las duraciones. La clase 0 % sale de los bloques de reposo, balanceada con 10 % y 40 %.
4. **Hurst** por canal. Port de `rangors.m` / `variograma.m` del grupo:
   - **HRS(p)**: R/S diádico con la regla de particiones de `rangors(v, p)`, p = 2…128 (Ecs. 7–9, 14–15).
   - **HO**: todas las particiones diádicas con bloque ≥ 8 (Ecs. 10–13).
   - **HV**: semivariograma γ(h), h = 1…N/2 (Ec. 16), al que luego se aplica R/S con p = 9, igual
     que en el script del grupo.
   - Si una regresión log-log tiene menos de 3 escalas válidas, el valor es NaN (`?` en ARFF).
5. **Weka.** Parámetros del artículo (§3.7.3):
   - **kNN**: IBk, k = 10, distancia de Chebyshev.
   - **MLP**: L = 0.3, M = 0.2, 500 épocas, 7 neuronas ocultas.
   - **SVM**: SMO con PolyKernel.
   - Opcionales con `cfg.weka.incluir_extra`: BN, NB, RF y RT.
   - Validación: CV estratificada de 5 folds con semilla 1. Métricas: % de aciertos, κ, MAE, RMSE,
     TPR, FPR, precisión, F, ROC y matriz de confusión.

### Diseño experimental (escalonado)

- **Datasets**: {KMI, FIS} × {3 clases, 2 clases} × {intra = una sesión, mes, global}.
- **Conjuntos de características**:
  - Tabla 2 del artículo: HRS64, HO, HV, HRS64+HO+HV.
  - Tabla 1 del artículo: HRSp+HO+HV para cada p.
- **Fase A — duración**: las 9 duraciones con los 16 canales (E9). La **mejor duración** se elige
  con un criterio fijado de antemano (`cfg.experimentos.criterio`): máxima precisión media
  intra-sujeto en KMI 3 clases con HRS64+HO+HV.
- **Fase B — canales**: 14 configuraciones en 512 muestras y en la mejor duración:

| ID | Canales | ID | Canales |
|---|---|---|---|
| E1 | F3, F4 | R1 | Cz |
| E2 | + Fz | R2 | Fz, Cz, Pz |
| E3 | + Fp1, Fp3 | R3 | Pz, P3, P4, P5, P6 |
| E4 | + P3, P4 | R4 | todos menos Cz (criterio del artículo) |
| E5 | + Pz | R5 | T7, T8, P8, O1, Oz (control negativo) |
| E6 | + Cz (= conf. 10 de la tesis hermana del laboratorio) | | |
| E7 | + P5, P6 | | |
| E8 | + T7, T8, P8 | | |
| E9 | los 16 | | |

- **Chequeo de fuga**: las 5 mejores combinaciones se reevalúan con CV **agrupada por bloque**.

## Etapa 6 — Ajuste de clasificadores (KMI, 3 clases)

`etapa6_ajuste_clasificadores.m`. Se aplica a las 13 sesiones intra-sujeto de KMI en 3 clases, con
512 muestras, configuración E9 y el conjunto HRS64+HO+HV. Todo se configura en `cfg.ajuste`.

1. **Diagnóstico de convergencia.**
   - MLP: precisión de entrenamiento vs CV por épocas y tasa de aprendizaje, con y sin parada temprana.
   - RF: precisión vs número de árboles.
   - Comparación de familias de Hurst.
2. **Búsqueda en rejilla** sobre kNN, SVM (polinomial y RBF), MLP y RF, cruzada con la **selección
   de características dentro de cada fold** (`AttributeSelectedClassifier`: ninguna, InfoGain top-8,
   16 o 24, CFS).
   - Se evalúa con CV 5-fold, semilla 1: la **misma CV** usada para reportar, por decisión del autor.
   - El criterio es **común a todas las sesiones**: un solo juego de hiperparámetros por
     clasificador, el de mayor precisión media.
3. **Ensamble**: `Vote`, con promedio de probabilidades de los 4 mejores modelos.
4. **Estimación final**: CV 5-fold **repetida 10 veces** para la línea base (parámetros del
   artículo), los modelos ajustados y el ensamble.
5. **Control de sesgo**: **CV anidada** (externa 5 × interna 5) para los 2 mejores clasificadores
   ajustados. Elegir los hiperparámetros sobre la misma CV es optimista; la cifra anidada es la
   estimación sin ese sesgo y conviene reportarla junto a la otra.
6. **Estadística**: binomial por sesión frente al azar (1/3); Wilcoxon frente al azar y del modelo
   ajustado frente a la línea base (Holm); prueba de permutación de etiquetas (100) en el dataset global.

Salidas en `resultados/ajuste/`:
- `tabla_ajuste_kmi3c.csv`: tabla para publicación.
- `hiperparametros_elegidos.csv`
- `busqueda_resumen.csv`
- `final_por_sesion.csv`
- `anidada_por_sesion.csv`
- `permutacion_global.csv`
- `figuras/`

## Etapa 7 — Calibración de caracterizadores

`etapa7_calibrar_caracterizadores.m`. HO y HRS(p) son la misma pendiente de log⟨R/S⟩ con distinto
rango de bloques. La etapa barre los 28 rangos válidos × 2 formas de promediar R/S y los pares de los
mejores, y elige el óptimo para estos datos. Incluye una CV anidada sobre la elección del rango.

Resultado: **Ro64a2 + Re128a2** (escalas cortas, bloques ≤ 64/128 muestras). HO y HRS64 del artículo
quedaron en los puestos 43 y 37 de 71.

## Etapa 8 — Reposo, asincronía de la KMI y sesgo de archivo

`analisis_confundido_reposo.m` mostró que, con las etiquetas de las etapas 2–7, el reposo se
clasificaba como el nivel de **su archivo** (72–80 %). El modelo reconocía la grabación, no el nivel.
Además, **la sesión imaginada fue asíncrona** (el paciente intentó el ejercicio cuantas veces quiso),
así que su `stims` **no** indica cuándo imaginó. Las etiquetas reposo/tarea de KMI de las etapas 2–7
no son válidas.

`etapa8_reposo_eventos.m` redefine reposo y tarea:

- **Tarea**:
  - **FIS** (síncrona): los bloques de `stims`.
  - **KMI** (asíncrona): intentos detectados por **ERD**. Es la caída de potencia mu/beta en Cz, Pz,
    P3 y P4 (`src/eventos/detectar_erd.m`). El detector se valida en FIS, donde los bloques se conocen.
- **Reposo**, dos definiciones, con ventanas de la misma longitud para todas las clases de un dataset:
  - **A**: primeros 5 s del registro (ventanas de 1 s).
  - **B**: [fin + 0.5 s, fin + 0.8 s] tras cada estímulo o intento (ventanas de 0.3 s).
- **Variante Δ**: se resta la media del reposo inicial del propio archivo para quitar su "huella".
- **Pruebas**:
  - `huella`: reposo del archivo de 10 % vs reposo del de 40 %. Lo ideal es ~50 %.
  - `tarea_rep`: tarea vs reposo, sin confundido de archivo.
  - `2c`: 10 vs 40 %.
  - `3c`: con un índice de **herencia del reposo**, el % de reposos clasificados como tarea que caen en
    el nivel de su archivo. Lo ideal es ~50 %.

Salidas en `resultados/eventos/`: `tabla_final_reposo.csv`, `validacion_erd_fis.csv`,
`eventos_por_archivo.csv` y `figuras/`.

## Etapas 9 y 10 — Algoritmos fractales adicionales (Niveles3) y calibración de Higuchi y DFA

- **Etapa 9** (`etapa9_fractales.m`, `corrida_fractales.m` → `resultados_asr/fractales/`):
  - Se portan los 8 algoritmos fractales de Niveles3 con sus parámetros:
    - RS con destendencia (`hurst_rs_tendencia.m`);
    - Higuchi (`higuchi_dimension.m`);
    - DFA1 (`dfa_exponente.m`);
    - pendiente/2 del variograma (`variograma_pendiente.m`);
    - PROM, HO, HRS64 y HV.
  - Con los parámetros por defecto, **solo Higuchi añade algo** a la base de Hurst (+0.9 puntos en
    3 clases). Juntar los 8 empeora el resultado (−4 puntos): son estimadores muy redundantes.
- **Etapa 10** (`etapa10_calibrar_hig_dfa.m`, `corrida_hig_dfa.m` → `resultados_asr/calib_hig_dfa/`):
  - Candidatos:
    - Higuchi en rangos de k (72);
    - DFA1 y DFA2 en rangos de escala (210).
  - Todos se calculan a partir de las curvas L(k) y F(n) de cada ventana (`curvas_ventana.m`,
    `caracteristicas_desde_curvas.m`).
  - Criterios de evaluación:
    - 3 clases dentro de la sesión, con CV anidada;
    - LOSO KMI 10 vs 40, con LOSO anidada (`loso_anidado.m`);
    - permutación.
  - Resultados:
    - En 3 clases ganan las **escalas muy cortas**: HIG k = 1–6 (56.0 %), DFA2 n = 4–20 (56.3 %) y
      el par HIG3a8 + DFAc4a24 (57.7 %; MLP 60.4 %). La CV anidada de la elección da ≈ 58.7 % (MLP).
      **Pero** la mejora viene de separar mejor 10/40 dentro de la sesión (88 %), mientras el reposo
      empeora (4 %) y su herencia sube de 79 % a 87 %. Es decir, se refuerza la huella del archivo.
    - En la generalización, el mejor individual (HIG k = 3–64: 62.6 %, p de permutación = 0.004) sale
      de elegir el mejor de 290 candidatos. La LOSO anidada, con la elección hecha dentro de cada
      fold, cae a 51 %. **No hay mejora robusta** sobre Re64a2 + Re128a4 (59.1 %).

## Etapa 11 — Conjuntos ganadores de Niveles3

`etapa11_niveles3_ganadores.m` (`corrida_niveles3_ganadores.m` → `resultados_asr/niveles3_ganadores/`)
reproduce en Weka los conjuntos de Niveles3 junto con sus clasificadores:

- **Conjuntos:** "Básicos per-channel" (RS + HIG + DFA + VAR), "Todos 7" y "Martínez".
- **Clasificadores:** LogReg con ridge = 1/(2C) y estandarización, SVM RBF con C = 12 y γ = 0.02, y
  RF con 300 árboles.

Se comparan con Hurst y con Higuchi + DFA calibrados.

| Evaluación | Conjunto | Exactitud |
|---|---|---|
| 10 vs 40 dentro de la sesión | N3 Básicos | 88 % |
| 10 vs 40 dentro de la sesión | Básicos calibrados | 92 % |
| Entre sesiones (LOSO) | N3 Básicos | 53–57 % |
| Entre sesiones (LOSO) | HIG3a8 + DFAc4a24 | 60.4 % (p = 0.058) |
| Entre sesiones (LOSO) | Re64a2 + Re128a4 | 59.1 % (p = 0.034) |

**Rendimiento:**
- 12 workers aprovechando los hilos SMT (`cfg.paralelo.workers`).
- Instancias de Weka en memoria en la LOSO (`weka_instancias.m`).
- Permutación secuencial de Besag–Clifford (`permutacion_secuencial.m`): se detiene con 20
  excedentes, como máximo 1000 permutaciones.

## Decisiones y supuestos (leer antes de interpretar resultados)

- **Montaje.** No es el Emotiv del artículo. Filas de `EEG_all`:
  `Fp1, Fp3, F3, Fz, F4, Cz, Pz, P5, P3, O1, Oz, T8, P8, P6, P4, T7`. No hay C3/C4.
  "Fp3" no es una posición estándar; probablemente sea **Fp2** (pendiente de confirmar con las notas
  de adquisición).
- **Sin EMG.** Los `.mat` solo traen 16 canales de EEG, así que la verificación de "contracción nula"
  del artículo no es reproducible. El QA usa la potencia de artefactos EEG como proxy.
- **Temporización.** `stims` no trae tiempos; se asume espaciado uniforme. La evidencia está en
  `qa_alineacion_estimulos.m`: en la ejecución física la potencia de artefactos sube en la mitad de
  tarea del ciclo.
- **Archivos excluidos.** Los 4 de Px.009 mes 1 y `px010_10_1` tienen formato *legado* (solo `B`,
  ~89 s). Su ciclo real de artefactos (~15–18 s) no coincide con el que implican sus `stims`, así que
  no se pueden etiquetar. `B` en los archivos normales es una copia de las últimas 2500 muestras de
  `EEG_all` y se ignora.
- **Nivel desde el nombre del archivo.** En Px006-M3 y Px009-M3 los cuatro archivos tienen el mismo
  valor en `stims` (40 y 10, respectivamente), y en Px008-M3-KMI están intercambiados respecto al
  nombre. Se usa el nombre, como en Niveles3; `stims` solo aporta la división reposo/tarea. Las
  discrepancias quedan como aviso en `registros.csv`. **Conviene confirmar Px008-M3-KMI** con la
  bitácora de adquisición.
- **R/S.** `rangors.m` escribe `mean(R / S)` con `/` matricial (cociente de mínimos cuadrados), no la
  media de R./S de la Ec. 13. Por defecto se usa la Ec. 13. Con
  `cfg.hurst.promedio_rs = 'rangors_original'` se reproduce el original (diferencia en H ≈ 0.003).
- **Diferencias con Niveles3**:
  - HV por R/S del semivariograma, no por la pendiente/2 del log-variograma.
  - Etiquetas reposo/tarea desde `stims`, no "todo el archivo = nivel".
  - ICA determinista por sesión con criterio de curtosis. El criterio anterior, "mayor varianza tras
    blanquear", es arbitrario porque ahí todas las varianzas valen ≈ 1.
  - Clasificación en Weka en lugar de scikit-learn.
- **Fuga temporal.** La CV estratificada por ventanas (como en el artículo) puede inflar la precisión
  porque ventanas del mismo bloque caen en entrenamiento y prueba. `reportes/chequeo_fuga.csv` lo
  cuantifica.
- **Comparaciones múltiples.** El "mejor de la rejilla" es optimista. Para la tesis, reporte la
  configuración prefijada (E9, 512 muestras) junto con la seleccionada.
- **Contribución de Cz.** El ICA limpia los 16 canales, así que las configuraciones sin Cz conservan
  una contribución mínima de Cz a través de la reproyección (solo en el subespacio de artefactos
  eliminado).

## Estructura

```
config_niveles.m  ejecutar_pipeline.m  etapa1..5_*.m  qa_alineacion_estimulos.m
src/io  src/preproc  src/segment  src/features  src/weka  src/reportes  src/util
pruebas/  (matlab.unittest; referencia/ = rangors.m y variograma.m originales)
Datos/    (entrada)      resultados/  (generado)
```
