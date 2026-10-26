function cfg = config_niveles()
%CONFIG_NIVELES  Configuración central del pipeline niveles4 (MATLAB + Weka).
%   cfg = config_niveles() devuelve una estructura con TODOS los parámetros
%   del pipeline. Ningún otro archivo contiene "números mágicos": cualquier
%   ajuste metodológico se hace aquí.
%
%   Referencia metodológica: Martinez-Peon et al. (2024) J. Neural Eng. 21
%   046024 — "Characterization and classification of kinesthetic motor
%   imagery levels". Migrado desde github.com/marg1875/Niveles3 (Python).

raiz = fileparts(mfilename('fullpath'));

%% ---------------------------------------------------------------- Rutas
cfg.rutas.raiz        = raiz;
cfg.rutas.datos       = fullfile(raiz, 'Datos');
cfg.rutas.resultados  = fullfile(raiz, 'resultados');
% Carpeta de resultados alternativa (p. ej. corrida con CAR): variable de
% entorno NIVELES_RESULTADOS con el nombre de la carpeta dentro de la raíz.
if ~isempty(getenv('NIVELES_RESULTADOS'))
    cfg.rutas.resultados = fullfile(raiz, getenv('NIVELES_RESULTADOS'));
end
cfg.rutas.weka_jar    = localizar_weka_jar();

%% ------------------------------------------------- Adquisición / montaje
cfg.fs = 250;                                   % g.Nautilus (Hz)
% Orden de las filas de EEG_all (Niveles3/config.py, ReporteAvances.tex).
% "Fp3" no es una posición estándar 10-20; probablemente sea Fp2.
cfg.canales.nombres = {'Fp1','Fp3','F3','Fz','F4','Cz','Pz','P5', ...
                       'P3','O1','Oz','T8','P8','P6','P4','T7'};

% Configuraciones de canales a comparar (E = escalonadas, R = por región).
E1 = {'F3','F4'};
E2 = [E1, {'Fz'}];
E3 = [E2, {'Fp1','Fp3'}];
E4 = [E3, {'P3','P4'}];
E5 = [E4, {'Pz'}];
E6 = [E5, {'Cz'}];                      % = conf. 10 de la tesis hermana
E7 = [E6, {'P5','P6'}];
E8 = [E7, {'T7','T8','P8'}];
E9 = cfg.canales.nombres;               % los 16
cfg.canales.configs = struct( ...
  'id',     {'E1','E2','E3','E4','E5','E6','E7','E8','E9', ...
             'R1','R2','R3','R4','R5'}, ...
  'canales',{E1, E2, E3, E4, E5, E6, E7, E8, E9, ...
             {'Cz'}, {'Fz','Cz','Pz'}, {'Pz','P3','P4','P5','P6'}, ...
             setdiff(E9, {'Cz'}, 'stable'), {'T7','T8','P8','O1','Oz'}}, ...
  'descripcion', {'Par DLPFC bilateral','+Fz (frontal/pre-SMA)', ...
             '+Fp1,Fp3 (prefrontal)','+P3,P4 (parietal)','+Pz (fronto-parietal sin motora)', ...
             '+Cz (= conf. 10 tesis hermana)','+P5,P6 (IPL lateral)', ...
             '+T7,T8,P8 (temporo-parietal)','Todos (16)', ...
             'Cz (referencia motora)','Linea media','Parietal / IPL', ...
             'Todos menos Cz (criterio del articulo)','Control negativo temporal/occipital'});
cfg.canales.config_referencia = 'E9';   % configuración usada en la fase A

%% ------------------------------------------------------ Inventario / QA
cfg.inventario.n_stims_esperado    = 80;
cfg.inventario.patron_ciclo        = [0 0 0 0 1 1 1 1];   % 0 = reposo, 1 = tarea
cfg.inventario.muestras_por_stim   = [470 510];           % rango válido N/K
cfg.inventario.excluir_formato_legado = false;  % los 5 legados entran al análisis 10 vs 40 (sin bloques de stims)

%% ------------------------------------------------------------- Filtrado
cfg.filtro.tipo          = 'butter';    % 'butter' (tesis/Niveles3) | 'fir' (artículo)
cfg.filtro.banda         = [0.5 45];    % Hz
cfg.filtro.orden_butter  = 4;
cfg.filtro.orden_fir     = 2*cfg.fs;    % orden par para fir1 pasabanda
cfg.filtro.notch_hz      = 60;
cfg.filtro.notch_q       = 30;

%% ------------------------------------- Referencia y rechazo de artefactos
% verificar_integridad.m mostró contaminación de modo común (en muchos
% archivos casi todos los canales correlacionan r > 0.98 tras filtrar 1-45 Hz)
% y muy distinta entre archivos -> referencia promedio común (CAR).
cfg.preproc.car               = true;
% Rechazo de ventanas: solo como análisis de sensibilidad (descarta ~50 %
% de las ventanas y no cambia la conclusión; ver comparar_variantes.m).
cfg.artefactos.activar        = false;
cfg.artefactos.umbral_z       = 6;      % |x - mediana| > 6 * 1.4826 * MAD en algún canal
cfg.artefactos.margen_s       = 0.25;   % se marca también +-0.25 s alrededor
cfg.artefactos.max_rechazo    = 0.5;    % aviso si un archivo pierde > 50 % de ventanas
cfg.artefactos.min_canales    = 1;      % canales que deben superar el umbral a la vez
% Electrodos puenteados: r > umbral entre canales tras CAR
cfg.puentes.metodo            = 'ninguno';   % 'ninguno' | 'eliminar' | 'interpolar'
cfg.puentes.umbral_r          = 0.99;
cfg.puentes.vecinos           = 3;
% Artifact Subspace Reconstruction (Mullen et al., 2015): repara ráfagas de
% alta varianza (picos) por archivo, tras puentes y antes del CAR, sin
% descartar ventanas. Ver asr_calibrar.m / asr_procesar.m.
cfg.asr.activar               = true;        % preprocesamiento final (comparar_asr.m)
cfg.asr.k                     = 30;          % corte en desviaciones (Chang et al. 2020: 20-30; 30 = más conservador)
cfg.asr.ventana_s             = 0.5;
cfg.asr.z_calib               = [-3.5 5.5];  % ventanas limpias para calibrar
cfg.asr.frac_canales          = 0.925;       % fracción de canales que deben estar en rango
cfg.asr.max_dims              = 0.66;        % máximo de componentes reparables

%% ------------------------------------------------- Blanqueo + ICA (Ecs. 2-6)
cfg.blanqueo.modo        = 'estandar';  % 'estandar' | 'p0_articulo'
cfg.blanqueo.tol_rango   = 1e-6;        % autovalores < tol*max se descartan
cfg.blanqueo.p0_pesos    = [0.8 0.5 0.2]; % mayores(2) / intermedios / menores(2)

cfg.ica.activar          = true;
cfg.ica.eta              = 0.05;        % tasa de aprendizaje inicial
cfg.ica.eta_max          = 1.0;         % tope del paso adaptativo
cfg.ica.max_iter         = 3000;
cfg.ica.tol              = 1e-3;        % paro: max|gradiente relativo|
cfg.ica.alfa             = 4;           % sub-gaussianas (kappa <= 0)
cfg.ica.beta             = 2;           % super-gaussianas (kappa > 0)
cfg.ica.max_reintentos   = 4;           % reintentos con eta/2 si diverge
cfg.ica.submuestreo      = 2;           % W se estima con 1 de cada n muestras
cfg.ica.curtosis_nsigma  = 2;           % artefacto: z robusto de log(kurt+3) > n
cfg.ica.max_removidos    = 2;

%% ---------------------------------------------------------- Segmentación
cfg.segmentacion.margen_inicio_s = 0.5;
cfg.segmentacion.margen_fin_s    = 0.25;
cfg.segmentacion.anclas_por_bloque = 3;
cfg.segmentacion.muestras_ranura   = 512;  % las anclas son centros de ranuras de 512
% Duraciones en MUESTRAS (512 = referencia del artículo; resto = segundos*fs)
cfg.segmentacion.duraciones = [512, round([2 1.5 1 0.5 0.4 0.3 0.2 0.1]*cfg.fs)];
cfg.segmentacion.duracion_referencia = 512;
cfg.segmentacion.balancear_reposo = true;

%% ---------------------------------------------------------------- Hurst
cfg.hurst.particiones_hrs  = [2 4 8 16 32 64 128];
cfg.hurst.p_hrs_principal  = 64;
cfg.hurst.ho_bloque_min    = 8;     % HO: todas las particiones con bloque >= 8
cfg.hurst.hv_particiones   = 9;     % HV = HRS(semivariograma, p=9) (script PM10)
cfg.hurst.min_escalas      = 3;     % mínimo de puntos para la regresión log-log
cfg.hurst.bloque_min_abs   = 2;     % ningún nivel con bloques < 2 muestras
% <R/S>_m: 'ecuacion13' (media de R./S, Ec. 13 del artículo) o
% 'rangors_original' (rangors.m usa '/' matricial: sum(R.*S)/sum(S.^2))
cfg.hurst.promedio_rs      = 'ecuacion13';

%% ------------------------------------- Algoritmos fractales de Niveles3 (etapa 9)
% Mismos parámetros que Niveles3/config.py (adaptativos a la longitud N).
cfg.fractales.rs_escala_min_s     = 0.128;   % RS: escala mínima = max(32, 0.128*fs)
cfg.fractales.rs_escala_max_frac  = 0.40;    %     escala máxima = 0.40*N
cfg.fractales.rs_n_escalas        = 20;      %     logspace
cfg.fractales.rs_min_bloques      = 5;
cfg.fractales.rs_min_escalas      = 5;
cfg.fractales.hig_kmax_frac       = 0.25;    % Higuchi: kmax = min(0.25*N, 30), >= 5
cfg.fractales.hig_kmax_tope       = 30;
cfg.fractales.hig_kmax_piso       = 5;
cfg.fractales.hig_min_puntos      = 5;
cfg.fractales.dfa_escala_min_s    = 0.064;   % DFA: escala mínima = max(16, 0.064*fs)
cfg.fractales.dfa_escala_max_frac = 0.33;    %      escala máxima = 0.33*N
cfg.fractales.dfa_n_escalas       = 15;
cfg.fractales.dfa_min_ventanas    = 3;
cfg.fractales.dfa_min_escalas     = 5;
cfg.fractales.var_lag_max_frac    = 0.20;    % Variograma: lag máx = max(0.2*N, 10), <= N/5
cfg.fractales.var_min_lags        = 10;
cfg.fractales.var_min_puntos      = 5;
cfg.fractales.familias = {'RS','HIG','DFA','VAR','PROM','HO','HRS64','HV'};

%% ------------------------------------- Calibración de Higuchi y DFA (etapa 10)
cfg.calib_hd.hig_kmin        = [1 2 3 4 6 8 12 16];
cfg.calib_hd.hig_kmax        = [6 8 10 12 16 20 24 30 40 48 64];
cfg.calib_hd.hig_min_puntos  = 5;          % valores de k en el rango
cfg.calib_hd.dfa_ns          = [4 5 6 8 10 12 16 20 24 32 40 48 64 80 96 128 170];
cfg.calib_hd.dfa_min_escalas = 4;          % escalas contiguas en el rango
cfg.calib_hd.min_puntos_validos = 4;       % puntos no NaN para la regresión
cfg.calib_hd.base_hurst      = {'Re16a2','Re64a4','Re64a2','Re128a4'};   % referencias
cfg.calib_hd.top_familia     = 6;          % finalistas por familia (pares y combinaciones)
cfg.calib_hd.n_permutaciones = 1000;

%% ------------------------------------- Laboratorio Higuchi/DFA (etapa 12)
cfg.lab.carpeta     = 'laboratorio';
cfg.lab.L           = [250 375 512 700 1000];   % 1, 1.5, 2, 2.8 (Niveles3) y 4 s
cfg.lab.hig_kmin    = [1 2 3 4 6 8 12 16];
cfg.lab.hig_kmax    = [6 8 10 12 16 20 24 32 48 64 96 128 192 250];   % se limita a L/4
cfg.lab.dfa_ns      = [4 5 6 8 10 12 16 20 24 32 40 48 64 80 96 128 170 250 333 500];
cfg.lab.dfa_ordenes = 1:3;                       % DFA1, DFA2, DFA3
cfg.lab.dfa_pasos   = [1 0.5];                   % ventanas internas sin / con 50 % de traslape
cfg.lab.top         = 5;                         % finalistas por familia y criterio
cfg.lab.bloque      = 25;                        % candidatos por trabajo (LOSO + excedente)

%% ------------------------------------- Complejidad: todas las familias (etapa 13)
cfg.complejidad.L          = [512 700 1000];
cfg.complejidad.hig_kmin   = [1 2 3 4 6 8 12 16];
cfg.complejidad.hig_kmax   = [6 8 10 12 16 20 24 32 48 64 96 128 192 250];
cfg.complejidad.ns         = [4 5 6 8 10 12 16 20 24 32 40 48 64 80 96 128 170 250 333];
cfg.complejidad.n_min      = [4 5 6 8 10 16 24 32];          % rejilla dispersa de rangos
cfg.complejidad.n_max      = [16 20 24 32 48 64 96 128 170 250 333];
cfg.complejidad.dfa_ordenes = 1:3;
cfg.complejidad.dfa_pasos  = [1 0.5];
cfg.complejidad.mf_ordenes = [1 3];
cfg.complejidad.mf_q       = [-5 -3 -1 1 3 5];
cfg.complejidad.mf_anchos  = [3 5];                          % h(-q) - h(q)
cfg.complejidad.ms         = [2 3 4 6 8 12 16 24 32 48 64 96 128 192 250];
cfg.complejidad.wavelets   = {'haar', 'db2', 'db4', 'sym4'};
cfg.complejidad.psd_f      = [1 2 3 4 6 8 10 13 16 20 25 30 35 40 45];
cfg.complejidad.box_s      = [4 8 16 32 64 128 256 512];
cfg.complejidad.se_m       = 1:3;  cfg.complejidad.se_r = [0.1 0.15 0.2 0.25 0.3];  cfg.complejidad.se_tau = [1 2 3 5 8];
cfg.complejidad.fe_r       = [0.1 0.2 0.3];                  cfg.complejidad.fe_tau = [1 2 4];
cfg.complejidad.pe_d       = 3:6;  cfg.complejidad.pe_tau = [1 2 3 4 6 8 12 16];
cfg.complejidad.lz_tau     = [1 2 4 8];
cfg.complejidad.mse_S      = [3 5 8];
cfg.complejidad.top_loso   = 30;                             % candidatos por familia y L que pasan a la LOSO
cfg.complejidad.adelante_pool = 12;                         % familias candidatas en la selección hacia adelante
cfg.complejidad.adelante_max  = 5;                          % familias como máximo en un conjunto
cfg.complejidad.n_perm     = 1000;
cfg.complejidad.max_finalistas = Inf;
cfg.complejidad.carpeta    = 'complejidad';

%% ------------------------------------- MFDFA h(+1) + Higuchi (etapa 14)
cfg.mfh.L        = 1000;                                     % ventanas de 4 s (curvas de la etapa 13)
cfg.mfh.familias = {'MFhp1', 'HIG', 'MFhp3', 'DFA3', 'WAV'};  % MF y HIG + controles
cfg.mfh.top      = 6;                                        % mejores de cada familia que se combinan
cfg.mfh.n_perm   = 1000;
cfg.mfh.carpeta  = 'mfdfa_higuchi';
cfg.mfh.dir_w    = fullfile('complejidad', 'cache');       % ventanas contiguas (curvas de la etapa 13)

%% ------------------------------------- Controles de revisión (etapa 15)
cfg.revision.carpeta = 'revision';
cfg.revision.dir_w   = fullfile('complejidad', 'cache');
cfg.revision.perm_todas = false;                   % permutación por sesión también para las configuraciones fractales
cfg.revision.dias_excluidos_orden = {};            % días sin hora real de grabación (control orden vs nivel)

%% ------------------------------------- Potencia de banda optimizada y combinada (etapa 16)
cfg.banda.carpeta      = fullfile('legados', 'banda_combinada');
cfg.banda.dir_w        = fullfile('legados', 'cache');       % ventanas de 4 s con los 56 registros
cfg.banda.fractales    = {'FEr10t4', 'FEr10t2', 'MFhp1o3n5a32', 'DFAto5a32', 'WVdb4j1a3', 'HIG12a96', 'MFhp3o3n5a24'};
cfg.banda.max_familias = 4;                                 % selección hacia adelante
cfg.banda.top_validar  = 5;                                 % mejores de B2 que pasan a validación completa
cfg.banda.n_perm       = 1000;                              % completas, sin paro temprano
cfg.banda.referencias  = {'ABSth+ABSmu+ABSb1+ABSb2', 'FEr10t2', 'MFhp1o3n5a32', 'DFAto5a32'};

%% ------------------------------------------------------------- Weka
cfg.weka.folds           = 5;
cfg.weka.semilla         = 1;
cfg.weka.svm_exponente   = 1.0;
cfg.weka.incluir_extra   = false;    % añade BN, NB, RF, RT (Tablas 1-3 completas)
cfg.weka.cv_modo         = 'estratificada';   % 'estratificada' | 'agrupada'
cfg.weka.top_fuga        = 5;        % combinaciones a re-evaluar con CV agrupada

%% ------------------------------------------------------- Experimentos
cfg.experimentos.condiciones = {'KMI','FIS'};
cfg.experimentos.esquemas    = {'3c','2c'};
cfg.experimentos.alcances    = {'intra','mes','global'};
cfg.experimentos.fase        = 'AB';      % 'A', 'B' o 'AB'
% Criterio prefijado de "mejor duración" (fase A)
cfg.experimentos.criterio.condicion = 'KMI';
cfg.experimentos.criterio.esquema   = '3c';
cfg.experimentos.criterio.alcance   = 'intra';
cfg.experimentos.criterio.conjunto  = 'HRS64+HO+HV';

%% ------------------------------- Etapa 6: ajuste de clasificadores (KMI 3c)
cfg.ajuste.condicion   = 'KMI';
cfg.ajuste.esquema     = '3c';
cfg.ajuste.duracion    = 512;
cfg.ajuste.config      = 'E9';
cfg.ajuste.conjunto    = 'HRS64+HO+HV';     % superconjunto; la selección elige dentro
cfg.ajuste.etiqueta    = '';                % sufijo de la carpeta de salida: resultados/ajuste<etiqueta>
cfg.ajuste.dir_arff    = '';                % '' = resultados/arff/d<duracion>; otra ruta = ARFF calibrados
cfg.ajuste.diag_conjuntos = {'HRS64','HO','HV','HRS64+HO','HO+HV','HRS64+HO+HV'};
cfg.ajuste.selecciones = {'ninguna','ig8','ig16','ig24','cfs'};
cfg.ajuste.rejilla.knn_K      = [1 3 5 7 10 15 20 30];
cfg.ajuste.rejilla.knn_dist   = {'euclid','manhattan','chebyshev'};
cfg.ajuste.rejilla.knn_peso   = {'ninguna','inv'};
cfg.ajuste.rejilla.svm_E      = [1 2 3];
cfg.ajuste.rejilla.svm_C_poly = [0.01 0.1 1 10 100];
cfg.ajuste.rejilla.svm_gamma  = [0.001 0.01 0.1 1];
cfg.ajuste.rejilla.svm_C_rbf  = [0.1 1 10 100];
cfg.ajuste.rejilla.svm_filtro = [0 1];       % 0 normaliza, 1 estandariza
cfg.ajuste.rejilla.mlp_H      = {'a','3','7','15','7,7'};
cfg.ajuste.rejilla.mlp_L      = [0.01 0.1 0.3];
cfg.ajuste.rejilla.mlp_N      = [250 500 1000];
cfg.ajuste.rejilla.mlp_V      = [0 20];      % % de validación para parada temprana
cfg.ajuste.rejilla.rf_arboles = [100 300 500];
cfg.ajuste.rejilla.rf_K       = [0 3 7 15];
cfg.ajuste.rejilla.rf_prof    = [0 10];
cfg.ajuste.semillas_final     = 1:10;        % CV 5-fold repetida 10 veces
cfg.ajuste.n_anidada          = 2;           % modelos con control de sesgo (CV anidada)
cfg.ajuste.anidada_max_config = 160;         % tope de la rejilla interna; si se excede se
                                             % submuestrea de forma equiespaciada (no
                                             % depende de los datos -> sin fuga)
cfg.ajuste.n_permutaciones    = 100;
cfg.ajuste.diag_mlp_N         = [10 50 100 250 500 1000 2000];
cfg.ajuste.diag_mlp_L         = [0.01 0.1 0.3];
cfg.ajuste.diag_rf_arboles    = [10 50 100 300 500 1000];

%% ------------------------- Etapa 7: calibración de caracterizadores (rangos R/S)
cfg.calibracion.top_pares       = 6;     % pares entre los 6 mejores candidatos individuales
cfg.calibracion.ajustar_despues = true;  % ejecutar la etapa 6 con el conjunto calibrado

%% ---------------------- Etapa 8: reposo y eventos (KMI asíncrona / FIS síncrona)
% La sesión imaginada (KMI) fue ASÍNCRONA: su stims no indica cuándo imaginó
% el paciente -> intentos por ERD. La física (FIS) es síncrona -> stims.
cfg.eventos.banda            = [8 30];             % mu/beta (Hz)
cfg.eventos.canales_erd      = {'Cz','Pz','P3','P4'};
cfg.eventos.factor_umbral    = 0.2;               % umbral = p10 + f (p90 - p10)
cfg.eventos.min_duracion_s   = 0.5;
cfg.eventos.separacion_min_s = 3.0;
cfg.eventos.reposo_inicial_s = 5.0;               % primeros 5 s = reposo seguro
cfg.eventos.reposo_desde_s   = 0.5;               % se evita el transitorio inicial
cfg.eventos.post_ini_s       = 0.5;               % reposo B: [fin+0.5, fin+0.8] s
cfg.eventos.L_A              = 256;               % definición A: ventanas de 1 s
cfg.eventos.L_B              = 75;                % definición B: ventanas de 0.3 s
cfg.eventos.anclas_fis       = 3;                 % ventanas de tarea por bloque FIS
cfg.eventos.min_por_clase    = 5;                 % mínimo para un dataset intra-sujeto
cfg.eventos.semillas         = 1:5;               % CV 5-fold repetida
% Caracterizadores (definidos por límites de bloque, válidos para cualquier L)
cfg.eventos.familias = struct( ...
    'nombre', {'Ro64a2', 'Re128a2', 'HO', 'HRS'}, ...
    'modo',   {'rangors_original', 'ecuacion13', 'ecuacion13', 'ecuacion13'}, ...
    'm_max',  {64, 128, Inf, Inf}, ...
    'm_min',  {2, 2, 8, 4});
cfg.eventos.conjuntos = {'Ro64a2+Re128a2', 'HO+HRS'};   % calibrado / referencia del artículo

%% ---------------------------------------------------- Ejecución
cfg.paralelo.activar   = true;
cfg.paralelo.workers   = 6;    % 8 núcleos: con 12-14 workers sostenidos el equipo se apagaba (corte de energía, Kernel-Power 41)
cfg.guardar_senales    = true;
cfg.modo_prueba        = false;   % true: solo Px006 mes 1, 2 duraciones, 2 configs
cfg.etapas             = [1 2 3 4 5 6];
cfg.version            = 'niveles4-1.0';
end

function jar = localizar_weka_jar()
% Busca weka.jar: variable de entorno WEKA_JAR o instalación en Program Files.
jar = getenv('WEKA_JAR');
if ~isempty(jar) && isfile(jar), return; end
cands = [dir('C:\Program Files\Weka-*\weka.jar'); dir('C:\Program Files (x86)\Weka-*\weka.jar')];
if isempty(cands)
    jar = '';
else
    jar = fullfile(cands(end).folder, cands(end).name);
end
end
