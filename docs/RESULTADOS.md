# niveles4 — Registro completo de resultados

Documento de resultados de todas las etapas (1–15) del pipeline MATLAB + Weka. Cada número
sale de los CSV y logs que se indican; las rutas son relativas a la raíz del repositorio.
La metodología y el código se describen en `README.md`. Las figuras generadas hasta ahora
son de trabajo (no de calidad de publicación) y se rehacerán en otra sesión.

**Abreviaturas**

| Abreviatura | Significado |
|---|---|
| KMI | Imaginería motora cinestésica (sesiones `*img*`), asíncrona |
| FIS | Ejecución física, síncrona por bloques de `stims` |
| 10 / 40 | Nivel de esfuerzo, % de la contracción voluntaria máxima |
| 3c | 3 clases: reposo / 10 % / 40 % |
| 2c | 2 clases: 10 % / 40 % |
| intra | CV de 5 folds dentro de una sesión (paciente × mes × condición) |
| LOSO | Dejar una sesión fuera (generalización a una sesión no vista) |
| LOPO | Dejar un paciente fuera |
| ig8 | Selección de 8 atributos por InfoGain dentro de cada fold |
| Huella | El clasificador reconoce el **archivo** en lugar del nivel |

---

## 0. Resumen ejecutivo

1. **Dentro de la sesión** los niveles 10/40 se separan con 70–92 % y 3 clases con 50–63 %.
   Gran parte de eso es **huella del archivo**:
   - el reposo hereda la clase de su archivo en un 72–94 % de los casos;
   - dos mitades del **mismo** archivo (mismo nivel) se distinguen al 62–67 %.
2. **Entre sesiones (LOSO)** el mejor resultado defendible es de **~60–67 %** en KMI 10 vs 40:
   - MFDFA h(+1) en ventanas de 4 s: 66.1 %, permutación p = 0.017;
   - MFDFA h(+1) + Higuchi: 66.8 %, p = 0.010;
   - la LOSO anidada estricta de toda la selección baja a **58–61 %** (no significativa).
3. **Por paciente (LOPO)** todo queda en 58–64 %, sin significancia (p ≥ 0.062; el p mínimo posible
   con 5 pacientes es 0.031).
4. **Por registro** (voto de las ventanas de un archivo), MFDFA h(+1) + Higuchi acierta el 85 % de
   los archivos y decide bien 11 de 13 sesiones (binomial p = 0.011).
5. **El orden de grabación no explica el resultado.** El registro KMI10, grabado después del FIS10,
   se predice *menos* "40", y el efecto del tiempo en el modelo mixto es ≈ 0.
6. **FIS no generaliza** entre sesiones con ninguna característica (43–58 %).
7. **El baseline estándar de potencia de banda θ/μ/β iguala a MFDFA** (LOSO 64.7 %, LOPO 63.4 %).
8. **Los resultados de Niveles3** (95 % en 2c, 85 % en 3c) se reproducen dentro de la sesión
   (88–92 %), pero no se sostienen entre sesiones (53–57 %).

---

## 1. Datos

Fuentes: `resultados/registros.csv`, `resultados/integridad/*.csv`, `resultados/qa/alineacion.csv`
y `resultados_car/puentes_por_archivo.csv`.

**Equipo y montaje**
- g.Nautilus de 16 canales a 250 Hz.
- Orden de canales: Fp1, Fp3 (≈ Fp2), F3, Fz, F4, Cz, Pz, P5, P3, O1, Oz, T8, P8, P6, P4, T7.
- Sin EMG y sin C3/C4.

**Registros**
- 5 pacientes con ECV (Px006–Px010), meses 1/3/6, KMI y FIS al 10 % y 40 %.
- 56 archivos: **51 válidos** (formato A: `EEG_all` + `stims` 1×80).
- 5 excluidos por formato legado (Px009 mes 1 ×4 y px010_10_1: solo `B`, ~89 s, sin temporización
  alineable).
- **26 sesiones** (25 completas con 10 y 40) → **13 sesiones KMI** y 12 FIS.

**Etiquetas de `stims`**
- `stims` es una secuencia de etiquetas `[0 0 0 0 L L L L] × 10` sin marcas de tiempo; se asume
  espaciado uniforme. El QA de alineación (`qa_alineacion_estimulos.m`) lo respalda: en FIS la
  potencia de artefactos sube en la mitad de tarea del ciclo.
- Inconsistencias:
  - En Px006-M3 todas las etiquetas valen 40 y en Px009-M3 todas valen 10; se usa el nivel del
    nombre de archivo.
  - En Px008-M3-KMI los nombres están intercambiados respecto a `stims`; se corrige con `stims` en el
    análisis temporal.
- `B` (en formato A) es una copia de las últimas 2500 muestras de `EEG_all`.

**Orden real de grabación** (marca "Created on" de los .mat):
- Casi siempre FIS10 → KMI10 (+~4 min) → FIS40 (+~9 min) → KMI40 (+~12–20 min).
- Excepción: **Px010-M1 KMI**, grabada en orden invertido (40 antes que 10).

**Calidad**
- Contaminación de modo común: en muchos archivos casi todos los pares de canales tienen r > 0.98
  después del filtro de 1–45 Hz.
- Canales duplicados en algunos archivos (bandera de `integridad_archivos.csv`).
- **Electrodos puenteados** (r > 0.99 tras CAR), concentrados en el **mes 6**:

  | Canales puenteados por archivo | Archivos |
  |---|---|
  | 0 | 17 |
  | 2–4 | 20 |
  | 5–8 | 9 |
  | 10–11 | 3 |

  Sesiones con más de 6 canales puenteados: Px006_M6_FIS, Px006_M6_KMI, Px007_M6_KMI,
  Px008_M6_FIS, Px008_M6_KMI y Px009_M6_FIS.

---

## 2. Preprocesamiento

**Pipeline final** (`resultados_asr/`):
- filtro DC + notch de 60 Hz + Butterworth de 0.5–45 Hz (orden 4, fase cero);
- ASR con k = 30;
- CAR;
- z-score;
- blanqueo PCA + ICA de gradiente natural con no linealidad flexible (una W por sesión);
- eliminación de hasta 2 fuentes por curtosis robusta;
- reproyección a electrodos.

**Estadísticas** (`resultados_var/asr_k30/pipeline.log`, `resultados_asr/pipeline.log`):

| Indicador | Valor |
|---|---|
| Fracción del registro reparada por ASR | 4.3 % (máx 30.7 %) |
| Pares de canales con r > 0.98 tras ASR + CAR + ICA, por archivo | 4.7 (máx 98) |
| Iteraciones medias del ICA | 1358 (máximo 3000) |
| **Sesiones que NO convergieron** (3000 iteraciones, tolerancia 1e-3) | **6 de 26**: Px006_M1_KMI, Px007_M1_KMI, Px007_M6_FIS, Px008_M1_KMI, Px008_M3_FIS, Px008_M6_KMI |
| Componentes eliminadas | 0–2 por sesión |

**Variantes de preprocesamiento comparadas** (la misma métrica, LOSO KMI 10 vs 40 con SVM ig8):

| Variante | LOSO KMI | p (permutación) | Ventanas usadas | Fuente |
|---|---|---|---|---|
| Sin CAR (etapas 1–8 originales) | ≈ 47 % | — | 100 % | `resultados/deriva` |
| CAR + rechazo de 6 MAD (P0-R1) | 58.4 % | 0.089 (100 perm.) | 51 % | `resultados_var/comparacion_variantes.csv` |
| Excluir sesiones con > 6 puentes + 6 MAD (P3-R1) | 59.5 % | 0.109 | 40 % | ídem |
| Interpolar puentes (mejor P1) | 58.1 % | — | 52 % | ídem |
| Eliminar puentes (mejor P2) | 58.4 % | — | 41 % | ídem |
| CAR sin rechazo (A0) | 56.4 % | 0.107 (1000 perm.) | 100 % | `resultados_var/comparacion_asr.csv`, `permutacion_final.csv` |
| ASR k = 20 | 55.8 % | — | 100 % | ídem |
| **ASR k = 30 (final)** | **59.1 %** | **0.034** (1000 perm.) | 100 % | ídem |

- **FIS en esas mismas variantes:** 43–54 %.
- **El rechazo de 6 MAD marca EEG normal, no artefactos:** después de ASR seguiría marcando ~50 % de
  las ventanas. Por eso es solo un análisis de sensibilidad.
- **Corregir los puentes** (interpolar o eliminar) no mejora frente a dejar los datos como están.
- **ASR con k = 30** se eligió entre k = 20 y k = 30; corrigiendo por esas dos opciones, p ≈ 0.07.

---

## 3. Etapas 1–5: réplica de Martinez-Peon et al. (sin CAR)

Fuente: `resultados/reportes/*.csv`. Características: Hurst HO, HRS(p) y HV por canal.
Clasificadores del artículo: kNN (k = 10, Chebyshev), MLP (L = 0.3, M = 0.2, N = 500, H = 7) y
SVM (PolyKernel). CV estratificada de 5 folds.

### 3.1 Fase A: duración de ventana (configuración E9 = 16 canales)

Exactitud por duración (máximo entre conjuntos y clasificadores / media). Las duraciones de 500 y 512
son casi idénticas; abajo se muestra 512.

| Alcance | Condición/esquema | 25 | 50 | 100 | 125 | 250 | 375 | 512 |
|---|---|---|---|---|---|---|---|---|
| intra | KMI 3c | 42.4/38.2 | 42.6/38.9 | 44.5/41.0 | 45.1/40.7 | 45.5/42.0 | 46.0/42.2 | **48.2/43.6** |
| intra | KMI 2c | 62.9/57.8 | 69.1/62.9 | 71.9/65.4 | 73.5/66.8 | 74.0/67.7 | 77.1/67.7 | 73.5/67.8 |
| intra | FIS 3c | 44.1/39.2 | 46.7/42.3 | 44.0/40.6 | 44.1/40.7 | 46.2/41.7 | 48.1/43.1 | **50.1/46.8** |
| intra | FIS 2c | 64.3/59.3 | 67.5/64.7 | 68.8/65.3 | 67.2/62.7 | 68.8/61.9 | 70.8/64.6 | **73.3/66.7** |
| global | KMI 3c | 39.7/35.7 | 39.0/36.3 | 37.7/35.9 | 40.2/35.6 | 38.7/35.9 | 40.3/36.5 | 39.1/36.4 |
| global | KMI 2c | 57.1/53.0 | 59.5/56.6 | 60.8/56.6 | 58.5/55.2 | 61.2/56.1 | 61.2/56.8 | **66.7/59.4** |

- **Las ventanas largas (2 s) son mejores** en casi todo.
- **Mezclar sesiones (global) reduce la exactitud casi al azar en 3 clases.**

### 3.2 Fase B: configuración de canales (512 muestras, intra)

Media de conjuntos × clasificadores, con el máximo entre paréntesis:

| Ranking | KMI 3c | KMI 2c | FIS 3c | FIS 2c |
|---|---|---|---|---|
| 1 | R4 43.8 (49.3) | **E9 67.8 (73.5)** | **E9 46.8 (50.1)** | **E9 66.7 (73.3)** |
| 2 | E9 43.6 (48.2) | E7 67.7 | E8 46.6 | R4 66.4 |
| 3 | E8 43.6 | R4 67.7 | R4 46.4 | E8 66.2 |
| último | R1 36.8 | R1 57.7 | R1 41.9 | R1 59.6 |

- **La progresión E1 → E9 es monótona:** más canales es mejor.
- **Friedman entre configuraciones** (13 sesiones): SVM p = 2.4e-5, kNN p = 0.034 y MLP p = 0.47
  (`faseB_friedman.csv`; el detalle por pares está en `faseB_wilcoxon_*.csv`).

### 3.3 Métodos del artículo: tablas 1–3 (512 muestras)

KMI, intra (13 sesiones):

| Conjunto | 3c kNN / MLP / SVM | 2c kNN / MLP / SVM |
|---|---|---|
| HO | 44.4 / 42.9 / **47.9** | 66.4 / 71.5 / **73.2** |
| HRS64 | 45.1 / 44.3 / **48.2** | 69.4 / 72.8 / **73.3** |
| HV | 39.6 / 42.9 / 41.6 | 55.4 / 69.6 / 66.3 |
| HRS64 + HO + HV | 41.3 / 43.8 / 45.1 | 57.9 / 71.3 / 71.9 |

- KMI global, 3c: 34–39 % (azar 33 %); 2c: 54–63 %.
- Las tablas de FIS y las de particiones (HRS2…HRS128) están en `tabla1_*`, `tabla2_*` y `tabla3_*`.
- **HV es el peor método** y combinar los tres no supera a HRS64 solo.

### 3.4 Fuga temporal (`chequeo_fuga.csv`)

La CV estratificada por ventanas (la del artículo) frente a la CV agrupada por bloque infla la
exactitud en promedio **+2.6 puntos** (2c +2.4 y 3c +2.8; 250 comparaciones, hasta +11.7).

### 3.5 Longitudinal por mes (`longitudinal_por_mes.csv`, intra)

| Condición | Mes 1 | Mes 3 | Mes 6 |
|---|---|---|---|
| KMI 2c (media / máx) | 69.2 / 77.9 | 67.9 / 75.3 | 66.2 / 74.2 |
| KMI 3c | 42.9 / 49.7 | 43.5 / 47.1 | 44.5 / 48.9 |
| FIS 2c | 56.6 / 64.4 | 64.7 / 71.3 | **76.7 / 86.2** |
| FIS 3c | 43.8 / 51.1 | 45.5 / 49.3 | 50.7 / 55.8 |

En FIS la exactitud sube con los meses. En KMI no hay progresión. Hay pocas sesiones por mes (3–5) y,
en el mes 6, electrodos puenteados.

---

## 4. Etapa 6: ajuste de clasificadores (KMI 3c, E9, 512 muestras)

Fuente: `resultados/ajuste*/tabla_ajuste_kmi3c.csv` y `hiperparametros_elegidos.csv`. Búsqueda
en rejilla con la misma CV, Vote, CV repetida y CV anidada.

| Conjunto | Mejor modelo | Exactitud (DE) | CV anidada | Sesiones > azar | Mejor modelo del artículo |
|---|---|---|---|---|---|
| HRS64 + HO + HV (base) | MLP ajustado | 49.4 % (10.6) | 47.8 % | 9/13 | RF 45.9 % / SVM 45.6 % |
| HRS64 + HO | MLP ajustado | 48.7 % | 47.8 % | 10/13 | SVM 47.4 % |
| **Calibrado** (Ro64a2 + Re128a2) | **MLP ajustado** | **52.5 %** (8.8) | **51.4 %** | **12/13** | SVM 50.9 % |

**Hiperparámetros elegidos** (conjunto calibrado):

| Clasificador | Hiperparámetros | Selección de atributos |
|---|---|---|
| SVM | Poly, C = 0.1, E = 3 | ig8 |
| MLP | H = 3, L = 0.01, N = 500, sin validación | ninguna |
| kNN | K = 30, euclidiana, sin pesos | ig8 |
| RF | K = 3, 100 árboles | ig16 |

**Permutación global** (MLP ajustado, 100 permutaciones): observado 39.3–40.5 % frente a una nula de
33.3 % ± 1.5; p = 0.0099 en los tres conjuntos.

**Diagnóstico de convergencia del MLP** (`diagnostico_convergencia.csv`): el entrenamiento sube con
las épocas (33 → 66 %), pero la CV se queda en ~33–39 %. Hay sobreajuste sin ganancia, de ahí H = 3 y
L = 0.01.

---

## 5. Etapa 7: calibración de los caracterizadores R/S

Fuente: `*/calibracion/calibracion_resumen.csv` y `calibracion_anidada.csv`. Los candidatos son la
pendiente de log⟨R/S⟩ en rangos contiguos de bloque diádico [m_max, m_min] (28 rangos × 2 formas
de promediar) y sus pares. Criterio: KMI 3c intra, media de 13 sesiones × 3 clasificadores.

| Preprocesamiento | Mejor conjunto | Exactitud | CV anidada (MLP) | HO del artículo (Re512a8) | HRS64 del artículo (Re512a4) |
|---|---|---|---|---|---|
| Sin CAR | Ro64a2 + Re128a2 | 51.5 % | 53.3 % | 47.3 % (puesto 43/71) | 47.7 % (37/71) |
| CAR + rechazo | Re64a2 + Re128a4 | 52.2 % | 54.0 % | 44.2 % (52/71) | 44.9 % (50/71) |
| ASR (final) | Re16a2 + Re64a4 | 53.5 % | 53.5 % | 43.0 % (51/71) | 43.4 % (49/71) |

- **Las escalas cortas** (bloques ≤ 64–128 muestras) superan a los HO/HRS del artículo, que usan todo
  el rango de escalas.
- **Los caracterizadores calibrados para 3c no son los que mejor generalizan.**
  - Re16a2 + Re64a4 da LOSO 55.5 %.
  - Re64a2 + Re128a4 da LOSO 59.1 %.

---

## 6. Etapa 8: reposo, asincronía de la KMI y sesgo de archivo

Fuente: `*/eventos/*.csv` y `*/analisis_confundido/resumen_reposo_por_archivo.csv`.

**El reposo hereda la clase de su archivo.** Ventanas de reposo de las etapas 2–7 (bloques 0 de `stims`):

| Escenario (sin CAR) | Reposo predicho como 10 en el archivo de 10 | Reposo predicho como 40 en el archivo de 40 | Herencia total |
|---|---|---|---|
| Calibrado / MLP ajustado | 71.9 % | 80.2 % | **78.7 %** (p ≈ 1e-271) |
| Artículo HRS64 + HO + HV / SVM | 56.9 % | 55.0 % | 75.6 % |

Con CAR la herencia fue 79–91 %; con ASR, 79–85 %.

**Detección de ERD** (para KMI asíncrona; validada en FIS, donde se conoce la tarea):

| Preprocesamiento | Eventos ERD dentro de bloques de tarea | Azar (fracción del tiempo en tarea) | Bloques detectados |
|---|---|---|---|
| Sin CAR | 57 % | 51.7 % | ~75 % |
| CAR | 59 % | 51.7 % | 85 % |
| ASR | 57.1 % | 51.7 % | 79.2 % |

El detector apenas supera al azar. En KMI se detectan ~18.7 intentos por archivo (rango 8–26).

**Redefinición del reposo** (A = primeros 5 s; B = 0.5–0.8 s tras cada estímulo o intento). Mejores
valores por prueba, sin CAR:

| Prueba | KMI A | KMI B | FIS A | FIS B |
|---|---|---|---|---|
| Huella (reposo del archivo 10 vs del archivo 40; ideal 50 %) | — | 59.8 % media (65.2 máx) | — | 64.6 % (72.5) |
| Tarea vs reposo del mismo archivo | 57.7 % | 56.6 % | 60.9 % | 50.2 % |
| 3c | 53.3 % (herencia 51 %) | 45.5 % (herencia 69 %) | 48.7 % | 41.9 % (herencia 76 %) |
| 2c | 75.6 % | 69.8 % | 72.4 % | 69.0 % |

- **Con ASR** la huella de KMI bajó a 57–59 %, la de FIS a 55–59 %, y la tarea frente al reposo fue
  de 52–65 %.
- **Con CAR + rechazo**, la definición A dio 3c de 85 %, pero con solo 5–6 sesiones: el rechazo
  eliminó el reposo de la mitad de las sesiones. Es un resultado no válido.
- **La variante Δ** (restar el reposo inicial del propio archivo) da hasta 88–93 %, pero **es
  inválida**: crea una huella artificial del archivo.

---

## 7. Deriva temporal, control de mitades y orden invertido

Fuente: `*/deriva/clasificacion_resumen.csv` y `prueba_orden_invertido.csv`.

**10 vs 40 dentro de la sesión frente al control de mitades** (mismo archivo, primera frente a
segunda mitad, ~1.3 min de separación):

| Preprocesamiento | KMI 10 vs 40 | KMI mitades | FIS 10 vs 40 | FIS mitades |
|---|---|---|---|---|
| Sin CAR | 80.7–84.1 % | 64.6–66.3 % | 76.5–80.6 % | 67.4–67.6 % |
| CAR | 83.3–86.5 % | 62.5–63.2 % | 74.7–80.8 % | 65.0–67.5 % |
| ASR | 80.9–85.2 % | 63.5–65.7 % | 74.1–80.0 % | 71.0–73.1 % |

- **El excedente del nivel sobre el control de mitades es de ~17–23 puntos.**
- **La señal no es estacionaria:** dos mitades del mismo nivel se distinguen al ~65 %.

**LOSO de las sesiones con orden normal** (z-score por sesión), con ASR y SVM ajustado: **57.1 %**
(DE 11.3; 12 sesiones).
- **Sesiones con exactitud muy baja** (posibles etiquetas o registros problemáticos):
  - Px007-M3-KMI: 8–13 %, sistemáticamente invertida;
  - Px008-M6-KMI: 12–31 %;
  - Px008-M1-KMI: 24–44 %.

**Px010-M1 (orden invertido):** se predice correctamente en todas las variantes.

| Preprocesamiento | Exactitud en Px010-M1 |
|---|---|
| Sin CAR | 57–79 % |
| CAR | 67–80 % (6/6 variantes) |
| ASR | 62–72 % |
| Ventanas de 4 s con MFDFA | 81–90 % |

**Análisis inválidos (descartados):**
- La destendencia lineal de la fatiga daba un 99 % artificial: la pendiente no extrapola.
- La normalización Δ, ya mencionada en la sección 6.

---

## 8. Etapa 9: algoritmos fractales de Niveles3 (2 s, ASR)

Fuente: `resultados_asr/fractales/*.csv`. Algoritmos: RS con destendencia, Higuchi, DFA1,
pendiente del variograma (VAR), PROM, HO, HRS64 y HV, con los parámetros por defecto de Niveles3.
Base: Re16a2 + Re64a4.

| Conjunto | 3c (media 3 clasificadores) | Solo / aporte sobre la base | LOSO KMI (SVM ig8) |
|---|---|---|---|
| Base | 53.3 % | — | 55.5 % |
| + HIG | **54.1 %** | 52.2 % / **+0.9** | 55.9 % |
| + DFA | 53.3 % | 47.5 % / 0.0 | 54.5 % |
| + VAR | 53.1 % | 48.9 % / −0.2 | 54.8 % |
| + HRS64 | 52.3 % | 43.4 % / −0.9 | 53.7 % |
| + HO | 52.0 % | 43.0 % / −1.2 | 53.4 % |
| + PROM | 52.0 % | 44.0 % / −1.3 | 55.3 % |
| + RS | 51.9 % | 42.1 % / −1.3 | 55.7 % |
| + HV | 50.1 % | 39.5 % / −3.1 | 55.5 % |
| Los 8 de Niveles3 | 47.9 % | — | 56.6 % (permutación p = 0.12) |
| Base + los 8 | 49.0 % | −4.3 | 56.1 % |

- **La selección hacia adelante solo agrega HIG** (CV anidada MLP: 55.7 %).
- **Matriz de confusión 3c (MLP):** reposo correcto 6.7 → 8.5 %, herencia 79 → 80 %.
- **Juntar algoritmos resta:** son estimadores redundantes del escalamiento.

---

## 9. Etapa 10: calibración de Higuchi y DFA (2 s, ASR)

Fuente: `resultados_asr/calib_hig_dfa/*.csv`. Candidatos: 72 de Higuchi (rango de k) y 210 de DFA
(rango de escalas, DFA1 y DFA2).

| Familia | Niveles3 por defecto: 3c / LOSO | Mejor 3c | LOSO de ese candidato | Mejor LOSO | 3c de ese candidato |
|---|---|---|---|---|---|
| Higuchi | 52.2 % / 55.1 % | HIG1a6 **56.0 %** (MLP 57.4) | 47.0 % | HIG3a64 62.6 % | 48.8 % |
| DFA | 47.5 % / 56.5 % | DFAc4a20 **56.3 %** (MLP 58.2) | 57.1 % | DFAc4a40 59.6 % | 54.8 % |

**Pares y combinaciones (3c):**
- HIG3a8 + DFAc4a24: **57.7 %** (MLP 60.4–60.7 %);
- Re16a2 + Re64a4 + DFAc4a24: 55.9 %.

**CV anidada 3c (MLP):** elegir entre Higuchi da 58.9 %, entre DFA 58.4 % y entre combinaciones 58.7 %.

**Confusión 3c (MLP):**

| Conjunto | Exactitud | Reposo | 10 % | 40 % | Herencia |
|---|---|---|---|---|---|
| Base | 55.0 % | 6.7 % | 79.5 % | 79.0 % | 79 % |
| HIG3a8 + DFAc4a24 | 60.4 % | 4.4 % | 88.2 % | 88.7 % | **86.7 %** |
| DFAc4a20 | 58.2 % | 5.9 % | 85.9 % | 82.8 % | 80 % |
| HIG1a6 | 57.4 % | 3.3 % | 85.6 % | 83.1 % | 85.4 % |

La mejora de 3c **refuerza la huella**: las escalas de 1–6 muestras son una firma del registro.

**Generalización (LOSO KMI) y permutaciones de 1000:**

| Conjunto | LOSO | p |
|---|---|---|
| HIG3a64 | 62.6 % | 0.004 (el mejor de 290: sesgo de selección) |
| Re64a2 + Re128a4 + DFAc4a20 + DFAc5a24 | 60.9 % | 0.045 |
| DFAc4a40 | 59.6 % | 0.081 |
| **LOSO anidada** (28 finalistas elegidos por 3c) | **51.1 %** | — |

---

## 10. Etapa 11: conjuntos ganadores de Niveles3 (2 s, ASR)

Fuente: `resultados_asr/niveles3_ganadores/*.csv`. Clasificadores de Niveles3 en Weka:
- LogReg (ridge 0.5 + estandarización);
- SVM RBF (C = 12, γ = 0.02);
- RF de 300 árboles.

Más los ajustados de este pipeline.

| Conjunto | Atributos | 3c intra | 2c intra | 2c por mes (1/3/6) | LOSO KMI (mejor modelo) |
|---|---|---|---|---|---|
| **HIG + DFA calibrados** | 32 | **60.7 %** (MLP) | 91.0 % (LogReg) | 97.1 / 87.0 / 90.0 | **60.4 %** (SVM ig8, p = 0.058) |
| Básicos calibrados (RS + HIG3a8 + DFAc4a24 + VAR) | 64 | 59.0 % | **92.3 %** (MLP) | 97.1 / 91.3 / 88.8 | 58.1 % (SVM RBF) |
| **N3 Básicos** (RS + HIG + DFA + VAR) | 64 | 55.6 % | 87.9 % (LogReg) | 89.6 / 89.7 / 84.2 | 57.0 % (ig8) / 53.1 % (LogReg) |
| Hurst 3c (Re16a2 + Re64a4) | 32 | 55.1 % | 84.7 % | 86.7 / 84.7 / 82.9 | 55.5 % |
| Hurst LOSO (Re64a2 + Re128a4) | 32 | 51.5 % | 84.9 % | 84.2 / 89.0 / 80.4 | **59.1 % (p = 0.034)** |
| N3 Todos 7 | 112 | 50.4 % | 83.6 % | 83.3 / 87.3 / 79.2 | 56.9 % |
| N3 Martínez (HO + HRS64 + HV) | 48 | 45.0 % | 70.8 % | 67.9 / 76.7 / 66.2 | 53.8 % |

- **Datos balanceados:** no se usó SMOTE.
- **Confusión de N3 Básicos (MLP):** reposo 6.9 %, herencia 80.8 %.
- **El conjunto ganador de Niveles3 se reproduce dentro de la sesión, pero no generaliza.** Su
  clasificador ganador (LogReg) es de los peores entre sesiones.

---

## 11. Etapa 12: laboratorio de Higuchi/DFA (longitud de ventana, orden y traslape)

Fuente: `resultados_asr/laboratorio/*.csv`. Longitudes de 250 a 1000 muestras (1–4 s), DFA de orden
1–3 con y sin traslape, Higuchi con k_max ≤ L/4; ~3700 candidatos.

| Ventana | Mejor 3c | Exactitud | Mejor LOSO | LOSO |
|---|---|---|---|---|
| 1 s | DFA3 con traslape, n = 5–32 | 54.3 % | HIG k = 12–16 | 61.8 % |
| 1.5 s | DFA2 con traslape, n = 4–24 | 56.5 % | HIG k = 6–20 | 61.7 % |
| 2 s | DFA2, n = 4–20 | 56.3 % | DFA3 con traslape, n = 5–32 | 63.8 % |
| **2.8 s** | **DFA3 con traslape, n = 5–32** | **57.2 %** (LOSO 63.2) | DFA3 con traslape, n = 5–24 | 64.2 % |
| **4 s** | DFA3, n = 5–32 | 56.5 % (LOSO 64.8) | **DFA3 con traslape, n = 5–32** | **64.8 %** |

**Anidado estricto** (elegir por 3c sin la sesión de prueba y tomar su LOSO):

| Conjunto de candidatos | LOSO | Sesiones > 50 % | p (Wilcoxon) |
|---|---|---|---|
| Todos | 63.3 % | 9/13 | 0.052 |
| Solo 4 s | **64.8 %** | 9/13 | **0.027** |
| Solo 2 s | 57.3 % | 9/13 | 0.24 |
| Regla de excedente | 48.3 % (falla) | — | — |

**Validación final** (rejilla de ventanas cada 4 s; permutación secuencial):

| Conjunto | LOSO | p | Observación |
|---|---|---|---|
| HIG1a64 + DFAt5a24 (4 s) | 67.4 % | 0.009 | Elegido mirando la LOSO: optimista |
| HIG3a64 (4 s) + DFAt5a32 (2.8 s) | 63.4 % | 0.027 | |
| HIG3a8 + DFAc4a24 (2 s) | 58.6 % | 0.14 | |
| Re64a2 + Re128a4 (2 s, en esta rejilla) | 56.5 % | 0.11 | |

La LOSO de una misma configuración varía ±3 puntos según dónde caigan las ventanas.

---

## 12. Etapa 13: 29 familias de caracterización (2, 2.8 y 4 s)

Fuente: `resultados_asr/complejidad/*.csv`. Familias implementadas en `src/complejidad/`:
- **Escalamiento:** Higuchi, DFA1–3 (con y sin traslape), DMA, MFDFA (h(q) con q = ±1, ±3, ±5 y
  ancho), varianza agregada, momentos absolutos, Hurst por wavelets (haar, db2, db4, sym4),
  exponente 1/f y box-counting.
- **Dimensiones directas:** Katz, Petrosian (2 variantes), Sevcik, Hjorth (movilidad y complejidad)
  y line length.
- **Entropías y complejidad:** SampEn (m, r, τ), FuzzyEn, PE, WPE, Lempel-Ziv (2 y 3 símbolos) y MSE.

Hubo ~5700 candidatos. Fases: 3c para todos → LOSO para los 30 mejores de cada familia → anidado →
selección hacia adelante anidada → validación con permutación.

### 12.1 Ranking por familia (mejor candidato por 3c y por LOSO)

Valores exactos de `ranking_familias.csv` (prefijo W = longitud de ventana en muestras: W512 = 2 s, W700 = 2.8 s, W1000 = 4 s).

| # | Familia | Mejor 3c | 3c | LOSO de ese candidato | Mejor LOSO | LOSO | 3c de ese candidato |
|---|---|---|---|---|---|---|---|
| 1 | FUZZYEN | W700FEr10t2 | 62.5 % | 61.0 % | W1000FEr10t4 | 65.7 % | 61.3 % |
| 2 | MFhp1 | W1000MFhp1o3n5a32 | 58.0 % | 66.1 % | W1000MFhp1o3n5a32 | 66.1 % | 58.0 % |
| 3 | DFA3 | W700DFAto5a32 | 57.4 % | 63.2 % | W1000DFAto5a32 | 64.8 % | 55.8 % |
| 4 | WPE | W512WPEd5t1 | 56.8 % | 53.0 % | W700WPEd3t8 | 58.8 % | 46.6 % |
| 5 | WAV | W512WVdb4j1a4 | 56.7 % | 58.4 % | W1000WVdb4j1a3 | 64.1 % | 55.0 % |
| 6 | MFhp3 | W1000MFhp3o3n5a32 | 56.5 % | 64.2 % | W1000MFhp3o3n5a24 | 64.4 % | 54.7 % |
| 7 | DFA2 | W512DFAc4a20 | 56.3 % | 57.1 % | W700DFAco6a16 | 61.8 % | 54.0 % |
| 8 | HIG | W512HIG1a6 | 55.8 % | 47.0 % | W1000HIG2a64 | 60.7 % | 49.8 % |
| 9 | PSD | W1000PSD8a45 | 55.7 % | 61.5 % | W1000PSD13a40 | 65.4 % | 52.2 % |
| 10 | PE | W1000PEd4t1 | 55.4 % | 50.1 % | W512PEd3t8 | 59.8 % | 45.5 % |
| 11 | PET | W512PETd | 55.4 % | 48.4 % | W512PETm | 55.8 % | 47.5 % |
| 12 | MFhp5 | W512MFhp5o3n6a32 | 54.8 % | 60.3 % | W700MFhp5o3n5a32 | 63.6 % | 53.3 % |
| 13 | BOX | W512BOX64a256 | 53.6 % | 41.7 % | W700BOX4a128 | 59.9 % | 48.5 % |
| 14 | DFA1 | W512DFAlo4a16 | 53.4 % | 47.6 % | W512DFAlo4a128 | 55.3 % | 50.2 % |
| 15 | DMA | W512DMA5a16 | 53.4 % | 48.3 % | W1000DMA6a250 | 57.8 % | 48.6 % |
| 16 | MFhm1 | W512MFhm1o3n6a64 | 53.0 % | 55.6 % | W1000MFhm1o3n8a48 | 61.0 % | 50.6 % |
| 17 | SAMPEN | W512SEm2r10t1 | 52.6 % | 50.2 % | W512SEm1r15t2 | 58.4 % | 50.3 % |
| 18 | REF | W512HIGn3 | 52.1 % | 55.0 % | W1000Re128a4 | 60.5 % | 47.5 % |
| 19 | HJORTH | W512HJmob | 51.1 % | 56.5 % | W512HJmob | 56.5 % | 51.1 % |
| 20 | MSE | W512MSE3 | 50.8 % | 53.8 % | W512MSE8 | 56.8 % | 46.4 % |
| 21 | SEV | W512SEV | 50.4 % | 56.6 % | W512SEV | 56.6 % | 50.4 % |
| 22 | LL | W512LL | 50.4 % | 55.9 % | W700LL | 56.5 % | 48.9 % |
| 23 | LZC | W512LZ3t1 | 49.7 % | 54.4 % | W700LZ2t4 | 59.4 % | 45.8 % |
| 24 | MFhm3 | W512MFhm3o3n6a96 | 49.0 % | 51.5 % | W512MFhm3o3n8a48 | 59.2 % | 47.5 % |
| 25 | MFhm5 | W1000MFhm5o1n16a333 | 48.6 % | 51.4 % | W512MFhm5o3n8a48 | 58.4 % | 45.6 % |
| 26 | KATZ | W512KATZ | 48.5 % | 53.2 % | W700KATZ | 56.6 % | 44.7 % |
| 27 | ABM | W1000ABM2a6 | 47.3 % | 55.5 % | W1000ABM2a32 | 59.9 % | 44.8 % |
| 28 | AGV | W512AGV2a6 | 46.8 % | 53.9 % | W512AGV3a24 | 59.8 % | 45.1 % |
| 29 | MFancho | W1000MFw3o1n10a96 | 42.6 % | 51.8 % | W1000MFw3o1n8a64 | 55.5 % | 41.2 % |

Las columnas "mejor LOSO" están infladas por elegir el máximo dentro de cada familia. Los valores
honestos son los del anidado.

### 12.2 Anidado estricto (elegir por 3c sin la sesión de prueba)

| Conjunto de candidatos | LOSO | Sesiones > 50 % | p (Wilcoxon) | Más elegido |
|---|---|---|---|---|
| MFhp1 | 66.1 % | 9/13 | 0.033 | W1000MFhp1o3n5a32 |
| MFhp3 | 64.1 % | 9/13 | 0.050 | W1000MFhp3o3n5a32 |
| W1000 | 63.5 % | 9/13 | 0.080 | W1000FEr10t2 |
| DFA3 | 63.2 % | 9/13 | 0.061 | W700DFAto5a32 |
| W512 | 62.3 % | 9/13 | 0.110 | W512FEr10t2 |
| **Todas las familias** | 62.2 % | 9/13 | 0.146 | W700FEr10t2 |
| FUZZYEN | 62.2 % | 9/13 | 0.146 | W700FEr10t2 |
| W700 | 61.0 % | 9/13 | 0.196 | W700FEr10t2 |
| MFhp5 | 60.6 % | 8/13 | 0.098 | W512MFhp5o3n6a32 |
| WAV | 58.4 % | 9/13 | 0.122 | W512WVdb4j1a4 |
| DFA2 | 58.4 % | 9/13 | 0.168 | W512DFAc4a20 |
| PSD | 58.2 % | 8/13 | 0.191 | W1000PSD8a45 |
| MFhm1 | 57.1 % | 10/13 | 0.127 | W512MFhm1o3n6a64 |
| SEV | 56.6 % | 8/13 | 0.089 | W512SEV |
| HJORTH | 56.5 % | 8/13 | 0.092 | W512HJmob |
| ABM | 55.5 % | 9/13 | 0.182 | W1000ABM2a6 |
| REF | 54.6 % | 7/13 | 0.414 | W512HIGn3 |
| LZC | 54.4 % | 8/13 | 0.191 | W512LZ3t1 |
| MSE | 54.2 % | 8/13 | 0.162 | W512MSE3 |
| AGV | 54.1 % | 8/13 | 0.250 | W512AGV2a6 |
| LL | 54.0 % | 9/13 | 0.144 | W512LL |
| KATZ | 53.2 % | 8/13 | 0.180 | W512KATZ |
| WPE | 53.0 % | 6/13 | 0.826 | W512WPEd5t1 |
| MFancho | 52.4 % | 8/13 | 0.243 | W1000MFw3o1n10a96 |
| MFhm3 | 51.6 % | 9/13 | 0.348 | W512MFhm3o3n6a96 |
| MFhm5 | 51.2 % | 7/13 | 0.457 | W1000MFhm5o1n16a333 |
| SAMPEN | 50.2 % | 8/13 | 0.722 | W512SEm2r10t1 |
| PE | 50.1 % | 7/13 | 0.839 | W1000PEd4t1 |
| PET | 49.4 % | 5/13 | 0.960 | W512PETd |
| HIG | 46.6 % | 4/13 | 0.331 | W512HIG1a6 |
| DMA | 45.7 % | 6/13 | 0.413 | W512DMA5a16 |
| DFA1 | 45.2 % | 5/13 | 0.376 | W512DFAlo4a16 |
| BOX | 41.7 % | 4/13 | 0.127 | W512BOX64a256 |

### 12.3 Selección hacia adelante entre familias (global y anidada)

| Ventana | Conjunto global | LOSO global | LOSO anidada |
|---|---|---|---|
| 2 s | PE + FuzzyEn + MFDFA h(+1) + WPE | 62.3 % | 62.0 % (10/13; p = 0.078) |
| 2.8 s | PE + FuzzyEn + Petrosian | 59.1 % | 59.1 % (p = 0.22) |
| 4 s | FuzzyEn + exponente 1/f | 63.2 % | 63.7 % (10/13; p = 0.106) |

Combinar familias no supera a MFDFA h(+1) sola.

### 12.4 Validación final (permutación secuencial, como máximo 1000)

| Conjunto | LOSO | p (permutación) | Sesiones > 50 % | Px010 | FIS |
|---|---|---|---|---|---|
| W1000MFhp1o3n5a32 | 66.1 % | 0.0170 (1000 perm.) | 9/13 | 81 % | 57.8 % |
| W1000FEr10t4 | 65.7 % | 0.0758 (264 perm.) | 9/13 | 78 % | 54.3 % |
| W1000MFhp1o3n5a20 | 65.6 % | 0.0170 (1000 perm.) | 9/13 | 90 % | 57.7 % |
| W1000DFAto5a32 | 64.8 % | 0.0269 (744 perm.) | 9/13 | 77 % | 55.6 % |
| W1000MFhp3o3n5a24 | 64.4 % | 0.0292 (684 perm.) | 9/13 | 87 % | 56.0 % |
| W1000MFhp3o3n5a32 | 64.2 % | 0.0238 (840 perm.) | 9/13 | 71 % | 54.8 % |
| W1000WVdb4j1a3 | 64.1 % | 0.0355 (564 perm.) | 9/13 | 77 % | 56.0 % |
| W700DFAto5a32 | 63.2 % | 0.0292 (684 perm.) | 9/13 | 79 % | 56.3 % |
| W1000FEr10t2+W1000PSD8a45 | 63.2 % | 0.0921 (228 perm.) | 10/13 | 81 % | 53.8 % |
| W512PEd6t1+W512FEr10t2+W512MFhp1o3n6a32+W512WPEd6t1 | 62.3 % | 0.0725 (276 perm.) | 10/13 | 74 % | 56.4 % |
| W700DFAco6a16 | 61.8 % | 0.0417 (480 perm.) | 9/13 | 69 % | 55.2 % |
| W700FEr10t2 | 61.0 % | 0.1389 (144 perm.) | 9/13 | 69 % | 55.7 % |
| W1000HIG2a64 | 60.7 % | 0.0309 (648 perm.) | 10/13 | 66 % | 53.2 % |
| W700PEd6t1+W700FEr10t2+W700PETd | 59.1 % | 0.1528 (144 perm.) | 9/13 | 69 % | 57.0 % |
| W512Re64a2+W512Re128a4 | 59.1 % | 0.0327 (612 perm.) | 10/13 | 68 % | 46.5 % |
| W700WPEd3t8 | 58.8 % | 0.0350 (600 perm.) | 8/13 | 67 % | 47.6 % |
| W512WVdb4j1a4 | 58.4 % | 0.1111 (180 perm.) | 9/13 | 70 % | 55.7 % |
| W512DFAc4a20 | 57.1 % | 0.1833 (120 perm.) | 8/13 | 76 % | 54.7 % |
| W512WPEd5t1 | 53.0 % | 0.3333 (60 perm.) | 6/13 | 64 % | 50.9 % |
| W512HIG1a6 | 47.0 % | 0.8333 (24 perm.) | 4/13 | 61 % | 52.1 % |

### 12.5 FuzzyEn: el detector de huella más potente

Confusión 3c (MLP):

| Conjunto | Exactitud | Reposo | 10 % | 40 % | Herencia |
|---|---|---|---|---|---|
| FuzzyEn W700 r = 0.10, τ = 2 | 63.5 % | **1.2 %** | 93.1 % | 96.2 % | **93.1 %** |
| FuzzyEn W512 | 64.3 % | 3.1 % | 94.9 % | 94.9 % | 94.1 % |
| MFDFA h(+1) W1000 | 59.4 % | 0.4 % | 90.4 % | 87.3 % | 89.2 % |
| DFA3 W700 | 59.0 % | 1.9 % | 91.5 % | 83.5 % | 85.8 % |
| Hurst Re64a2 + Re128a4 | 51.5 % | 4.1 % | 74.1 % | 76.4 % | 81.3 % |

---

## 13. Etapa 14: MFDFA h(+1) + Higuchi (4 s)

Fuente: `resultados_asr/mfdfa_higuchi/*.csv`. Se evaluaron 562 variantes (MFDFA h(+1), Higuchi y,
como controles, MFDFA h(+3), DFA3 y wavelets) y 66 combinaciones.

**Mejores combinaciones** (SVM ig8, 1000 permutaciones):

| Configuración | LOSO | p | Sesiones > 50 % | Px010 | FIS |
|---|---|---|---|---|---|
| **MFDFA h(+1) o3 n6–24 + HIG k12–32** | **66.8 %** | **0.010** | **11/13** | 87 % | 57.4 % |
| MFDFA h(+1) o3 n5–24 + HIG k4–48 | 66.8 % | 0.011 | 10/13 | 83 % | 56.0 % |
| MFDFA h(+1) o3 n6–24 + HIG k3–64 | 66.7 % | 0.010 | 10/13 | 86 % | 56.5 % |
| MFDFA h(+1) o3 n5–24 + HIG k3–64 | 66.5 % | 0.012 | 9/13 | 84 % | 56.1 % |
| MFDFA h(+1) + HIG + MFDFA h(+3) | 66.4 % | 0.012 | 11/13 | 86 % | 54.4 % |
| MFDFA h(+1) o3 n5–32 sola | 66.1 % | 0.017 | 9/13 | 81 % | 57.8 % |
| HIG k3–64 sola | 62.9 % | 0.007 | 11/13 | 66 % | 51.9 % |

- **Hay una meseta:** unas 10 combinaciones dan entre 65.8 % y 66.8 %.
- **El clasificador importa:** el mejor par con SVM ig16 da 61–64 % y sin selección 54–57 %.

**LOSO anidada estricta de todo el procedimiento** (se repite la elección completa sin cada sesión):

| Regla | LOSO | Sesiones > 50 % | p (Wilcoxon) | Configuración más elegida |
|---|---|---|---|---|
| Solo MFDFA | **60.6 %** | 7/13 | 0.105 | o3 n5–24 (4/13) |
| MFDFA + Higuchi | 58.4 % | 8/13 | 0.151 | variable (1/13) |
| Solo Higuchi | 55.4 % | 9/13 | 0.152 | k3–64 (3/13) |

---

## 14. Etapa 15: controles de revisión

Fuente: `resultados_asr/revision/*.csv`. Rejilla de ventanas de 4 s (Hurst en 2 s), KMI. Intervalos
de confianza del 95 % por bootstrap de pacientes. El LOPO usa permutación exacta por paciente
(32 asignaciones).

### 14.1 LOSO, LOPO y voto por registro

| Configuración | LOSO [IC 95 %] | LOPO [IC 95 %] | p por paciente | Voto por archivo (26) | Decisión por sesión (13) |
|---|---|---|---|---|---|
| **MFDFA h(+1) + Higuchi** | 66.8 % [56–78] | 61.6 % [50–72] | 0.19 | **85 %** (p < 0.001) | **85 %, 11/13** (p = 0.011) |
| MFDFA h(+1) n5–32 | 66.1 % [49–82] | 61.7 % [48–73] | 0.19 | 81 % (p = 0.001) | 69 % (p = 0.13) |
| MFDFA h(+1) n5–24 | 64.6 % [49–80] | 62.1 % [48–75] | 0.16 | 81 % (p = 0.001) | 69 % (p = 0.13) |
| DFA3 con traslape | 64.8 % [47–81] | 60.7 % [46–74] | 0.16 | 73 % (p = 0.014) | 69 % (p = 0.13) |
| Hurst Re64a2 + Re128a4 | 59.1 % [57–62] | 58.1 % [55–62] | 0.062 | 77 % (p = 0.005) | 77 % (p = 0.046) |
| FuzzyEn | 63.3 % [53–75] | 64.4 % [53–77] | 0.062 | 73 % | 69 % |
| **Potencia θ/μ/β (SVM ig8)** | **64.7 %** [53–79] (perm. p = 0.096) | **63.4 %** [51–79] | 0.062 | 69 % | 77 % (p = 0.046) |
| Potencia θ/μ/β (LogReg) | 42.1 % | 43.6 % | 0.75 | 35 % | 38 % |
| CSP 8–30 Hz + LogReg | 58.1 % (perm. p = 0.28) | 54.6 % | 0.31 | 58 % | 54 % |

Los archivos y sesiones de un mismo paciente no son independientes, así que la binomial del voto es
optimista.

### 14.2 Orden temporal frente a nivel

Se entrena con la KMI de los otros días y se predicen los 4 registros de cada día (12 días). La tabla
da la fracción de ventanas predichas como "40" y un modelo mixto
`frac40 ~ nivel40 + minuto + FIS + (1 | día)`.

| Configuración | KMI10 − FIS10 (mismo nivel, ~4 min después) | KMI40 − KMI10 | FIS40 − FIS10 | Efecto del nivel (p) | Efecto del tiempo (p) |
|---|---|---|---|---|---|
| MFDFA h(+1) + Higuchi | **−0.16** | +0.24 | +0.01 | +0.12 (**0.018**) | ≈ 0 (0.40) |
| MFDFA h(+1) n5–24 | −0.05 | +0.08 | +0.03 | +0.05 (0.35) | ≈ 0 (0.32) |
| MFDFA h(+1) n5–32 | −0.07 | +0.08 | 0.00 | +0.03 (0.59) | ≈ 0 (0.26) |
| DFA3 | −0.05 | +0.08 | +0.01 | +0.03 (0.50) | ≈ 0 (0.11) |
| Hurst | −0.07 | +0.14 | +0.06 | +0.09 (0.015) | ≈ 0 (0.31) |
| FuzzyEn | −0.06 | +0.25 | **+0.21** | +0.21 (0.007) | ≈ 0 (0.055) |
| Potencia de banda | −0.04 | +0.24 | **+0.19** | +0.20 (0.021) | ≈ 0 (0.15) |

- **En todas las configuraciones KMI10 se ve menos "40" que FIS10**, que se grabó antes. Es lo
  contrario de un clasificador de "más tarde en el día".
- **El nivel se transfiere a la física solo con potencia de banda y FuzzyEn.**

---

## 15. Conclusiones para la publicación

**Defendible**
- Preprocesamiento estándar (ASR + CAR + ICA, sin descartar datos).
- Información del nivel de esfuerzo imaginado (10 vs 40 % de la MVC) que generaliza a sesiones no
  vistas: ~60–67 % por ventana y ~85 % por registro.
- La información está en el escalamiento de 20–130 ms con ventanas largas (MFDFA h(+1), DFA3,
  wavelets).
- Controles de huella del archivo (herencia del reposo, control de mitades) como contribución
  metodológica.
- El control de orden temporal no apoya la explicación por fatiga u orden.

**No defendible**
- La clasificación de 3 clases y las cifras dentro de la sesión como rendimiento real.
- La superioridad de los fractales sobre la potencia de banda.
- Significancia a nivel de paciente.
- Los resultados de Niveles3 como generalizables.

**Limitaciones**
- 5 pacientes y 13 sesiones KMI.
- Nivel confundido con el archivo y orden casi fijo (solo una sesión invertida).
- KMI asíncrona con etiquetas por archivo.
- Electrodos puenteados en el mes 6.
- 6 sesiones sin convergencia del ICA en 3000 iteraciones.
- Miles de configuraciones exploradas: reportar el anidado.

**Siguientes pasos recomendados**
1. Datos nuevos con el orden de 10/40 contrabalanceado.
2. Congelar la configuración (MFDFA h(+1) orden 3, n = 5–24/32, ventanas de 4 s, SVM ig8) y
   validarla prerregistrada.
3. Comparar y combinar MFDFA con potencia de banda en el voto por registro.
4. Revisar las sesiones atípicas (Px007-M3-KMI con LOSO de 8–13 %, invertida de forma sistemática;
   Px008-M6-KMI) y la no convergencia del ICA.
5. Rehacer las figuras con calidad de publicación.

---

## 16. Mapa de archivos de resultados

| Carpeta | Contenido |
|---|---|
| `resultados/` | Etapas 1–8 sin CAR (réplica del artículo) |
| `resultados/reportes/` | Fases A/B, tablas 1–3, fuga, longitudinal, Friedman/Wilcoxon |
| `resultados/ajuste*/` | Etapa 6 (base, HRS64 + HO, calibrado) |
| `resultados/integridad/`, `resultados/qa/` | Integridad de archivos y canales; alineación de `stims` |
| `resultados_car/` | Etapas 1–8 con CAR + rechazo; `puentes_por_archivo.csv` |
| `resultados_var/` | Variantes de puentes × rechazo y ASR k = 20/30; permutaciones |
| `resultados_asr/` | Pipeline final (ASR k = 30) y etapas 7–15 |
| `resultados_asr/fractales/` | Etapa 9 |
| `resultados_asr/calib_hig_dfa/` | Etapa 10 |
| `resultados_asr/niveles3_ganadores/` | Etapa 11 |
| `resultados_asr/laboratorio/` | Etapa 12 |
| `resultados_asr/complejidad/` | Etapa 13 (`ranking_familias.csv`, `ranking_candidatos.csv`, `anidado.csv`, `adelante.csv`, `validacion_final.csv`) |
| `resultados_asr/mfdfa_higuchi/` | Etapa 14 |
| `resultados_asr/revision/` | Etapa 15 |

Los archivos pesados (`.mat`, `.arff`, cachés, características por ventana) y los datos de
pacientes no se versionan; se regeneran con los scripts `corrida_*.m` (ver README).
