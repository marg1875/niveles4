# Experimento: controles sanos con el pipeline de los pacientes

> Rama `pacientes-control`. Experimento exploratorio (31-oct-2026); **no forma parte del artículo**.
> Puede servir para la tesis.

## 1. Participantes

Cinco controles sanos registrados con el mismo protocolo de niveles4 (FIS y KMI al 10 % y 40 % de la
MVC), por el mismo experimentador y con el mismo g.Nautilus de 16 canales a 250 Hz. Hay una sesión por
persona. Se eligieron de una lista de 12 voluntarios para comparar el **grupo** de pacientes contra el
**grupo** de controles (no por parejas):

| | Pacientes (Px006–Px010) | Controles (C01–C05) |
|---|---|---|
| Sexo | 4 H / 1 M | 4 H / 1 M |
| Edad, media ± DE | 50.4 ± 13.0 | 49.6 ± 14.8 |
| Mediana (rango) | 53 (32–64) | 53 (33–68) |
| Mano de la prensión | 3 derecha / 2 izquierda (la menos afectada) | 5 derecha |

- Edad: sumas de rangos idénticas (Mann-Whitney p = 1).
- La lateralidad previa al EVC de Px007 y Px010 (los que usaron la mano izquierda) no se conoce.
- Los datos están en `Datos_controles/` (no se suben). Las claves C01–C05 son anónimas; la
  correspondencia está solo en `Datos_controles/clave_controles.csv`.

## 2. Datos

- Mismo formato que los pacientes: `EEG_all` 16 × ~39 200, `stims` 1 × 80 con patrón
  `[0 0 0 0 L L L L] × 10`, N/K = 488–493. `listar_registros` acepta los 20 archivos.
- Orden de grabación en los 5 controles: FIS10 → KMI10 → FIS40 → KMI40.
- Cada zip traía además un archivo largo (~644 s, solo `B`, 320 dígitos; se responde solo al "3")
  de otra tarea de atención; no se usó.
- **Calidad:** canales sin señal (≈ 0 µV en 1–45 Hz):

  | Canal | Controles (de 20 archivos) | Pacientes (de 51) |
  |---|---|---|
  | Oz | 17 | 1 |
  | PO8 | 13 | 1 |
  | T8 | 12 (todo C01, C02 y C03) | 0 |
  | Cz | 6 (todo C05, C01 al 40 %) | 4 |
  | P3 / Pz | 4 / 4 | 0 / 8 |
  | F3, Fz, C3, C4 | 0 | ≤ 3 |

  Además, Fp1 de C01 tiene 5–12 mV (sin contacto) y en C02–C04 domina el modo común.

## 3. Método

`corrida_controles.m` → `resultados_controles/`. Mismas funciones y parámetros que los pacientes,
con los 16 canales:

1. Etapa 1: filtro, ASR k = 30, CAR e ICA por sesión.
2. Tres clases dentro de la sesión (reposo, 10 %, 40 %): HRS64 + HO + HV, 512 muestras, SVM del
   artículo, CV de 5 particiones con 10 semillas (escenario de `fig_confusion`).
3. Ventanas contiguas de 4 s y las 9 configuraciones de la tabla final de los pacientes
   (`validacion_final.csv`) con el SVM ajustado (Poly E = 3, C = 0.1, InfoGain-8): dentro de la
   sesión (5 particiones × 10 semillas), LOSO, voto por registro y decisión por sesión.
   - Con una sesión por persona, la LOSO de los controles es también la LOPO: se compara con la
     **LOPO** de los pacientes.
   - p: permutación exacta de las 2⁵ = 32 combinaciones de sesiones volteadas (mínimo 0.031).
4. Las cifras de los pacientes se recalcularon con el mismo código y coinciden con las publicadas.

Preprocesamiento: el ICA no convergió en 4 de 10 sesiones (C01 KMI y FIS, C03 FIS, C05 FIS);
ASR reparó 0–32 % de cada registro.

## 4. Resultados (KMI 10 % vs 40 %)

| Configuración | Intra pac. | Intra ctl. | LOSO pac. | LOPO pac. | Entre personas ctl. | p ctl. | Voto pac. / ctl. | Sesión pac. / ctl. |
|---|---|---|---|---|---|---|---|---|
| FuzzyEn τ4 + τ2 + Higuchi | 92.7 | 82.8 | 72.8 | 68.6 | 41.8 | 0.75 | 85.7 / 50 | 85.7 / 20 |
| FuzzyEn + potencia γ, β1 | 96.1 | 83.2 | 71.3 | 67.1 | 36.2 | 0.88 | 78.6 / 40 | 78.6 / 20 |
| Potencia θ + μ + β1 + β2 | 90.8 | 85.5 | 67.3 | 66.0 | 36.7 | 0.88 | 82.1 / 40 | 78.6 / 20 |
| FuzzyEn τ4 | 91.7 | 77.7 | 66.7 | 66.4 | 43.3 | 0.75 | 75.0 / 40 | 71.4 / 40 |
| FuzzyEn τ2 | 91.4 | 77.8 | 66.3 | 65.4 | 43.1 | 0.88 | 78.6 / 50 | 78.6 / 0 |
| Potencia, 6 bandas | 92.9 | 84.5 | 65.9 | 64.1 | 37.4 | 0.84 | 75.0 / 40 | 78.6 / 0 |
| Potencia μ | 89.9 | 82.9 | 63.9 | 60.2 | 48.5 | 0.59 | 75.0 / 40 | 71.4 / 60 |
| DFA3 | 88.7 | 80.2 | 62.2 | 55.9 | 41.8 | 0.84 | 64.3 / 30 | 64.3 / 20 |
| MFDFA h(+1) | 89.0 | 82.9 | 61.1 | 55.6 | 42.1 | 0.81 | 71.4 / 30 | 64.3 / 40 |

Voto: % de archivos (28 en pacientes, 10 en controles). Sesión: % de sesiones (14 / 5) en las que el
registro de 40 % recibe más ventanas predichas como 40 %.

Mejor configuración, por control:

| Control | Dentro de la sesión | Entre personas |
|---|---|---|
| C01 | 52.8 | 51.3 |
| C02 | 91.5 | 35.9 |
| C03 | 97.9 | 33.3 |
| C04 | 97.4 | 41.0 |
| C05 | 74.5 | 47.4 |

### Matrices de confusión

10 vs 40, mejor configuración (% de cada clase real, media entre sesiones):

| | Pac. intra | Pac. LOSO | Ctl. intra | Ctl. entre personas |
|---|---|---|---|---|
| Real 10 → pred. 10 / 40 | 88.5 / 11.5 | 66.1 / 33.9 | 85.4 / 14.6 | 32.3 / 67.7 |
| Real 40 → pred. 10 / 40 | 3.2 / 96.8 | 20.4 / 79.6 | 19.7 / 80.3 | 48.7 / 51.3 |

Tres clases dentro de la sesión (HRS64 + HO + HV, SVM del artículo): pacientes 44.4 ± 8.1 %
(13 sesiones), controles 45.9 ± 10.7 % (5 sesiones); azar 33.3 %.

| Real → pred. (%) | Reposo | 10 % | 40 % |
|---|---|---|---|
| Pacientes, reposo | 33.9 | 35.1 | 31.0 |
| Pacientes, 10 % | 35.1 | 49.8 | 15.2 |
| Pacientes, 40 % | 30.4 | 19.9 | 49.6 |
| Controles, reposo | 28.6 | 37.3 | 34.1 |
| Controles, 10 % | 26.9 | 58.3 | 14.9 |
| Controles, 40 % | 30.2 | 19.1 | 50.7 |

## 5. Lectura

- **Dentro de la sesión** los controles separan los niveles algo menos que los pacientes
  (78–86 % frente a 89–96 %). La diferencia la produce sobre todo C01 (53 %); sin C01 la mejor
  configuración da 90.3 %, como en los pacientes. El número de ventanas por sesión es el mismo
  (~78 en controles, ~75 en pacientes).
- **Entre personas** los controles no generalizan: 36–49 %, por debajo del azar y sin significancia.
  La diferencia 10 → 40 va en sentidos opuestos según la persona (C02 se predice invertido). Los
  pacientes sí generalizan entre personas (LOPO 56–69 %).
- **Tres clases**: igual en ambos grupos (~45 %).
- **Posibles causas** (no comprobadas):
  - Cantidad de datos: cada modelo de controles se entrena con 4 sesiones (~310 ventanas); la LOPO de
    los pacientes, con ~11 sesiones (~840 ventanas) de 4 personas. Prueba pendiente: repetir la LOPO
    de los pacientes con una sola sesión por paciente (162 combinaciones).
  - Canales sin señal distintos en cada control.
  - El 10 % siempre se grabó antes que el 40 %.

## 6. Archivos

- `corrida_controles.m`: la corrida completa.
- `resultados_controles/legados/banda_combinada/tabla_pacientes_controles.csv`, `confusion_2c.csv`,
  `por_sesion.csv`.
- `resultados_controles/analisis_confundido/confusion_3c.csv`.
- `src/reportes/predicciones_sesion.m` (antes función local de `analisis_confundido_reposo.m`) y
  `src/reportes/intra_por_sesion.m` (devuelve también la matriz de confusión).
