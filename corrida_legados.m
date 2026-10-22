%% CORRIDA_LEGADOS  Resultados del artículo con los 5 registros legados incluidos.
%  Los registros de formato legado (Px009 mes 1: FIS10, FIS40, KMI10, KMI40;
%  px010_10_1) se habían excluido por no tener temporización de estímulos.
%  El análisis principal (10 vs 40 con la etiqueta del archivo) no la
%  necesita, así que aquí se incluyen: 14 sesiones KMI y 14 FIS (antes 13 y 12).
%  1) Ventanas contiguas de 2 s y 4 s: se reutilizan las de la etapa 13 y
%     se agregan las sesiones nuevas o recompletadas (Px009_M1_KMI,
%     Px009_M1_FIS, Px010_M1_FIS), con las mismas curvas (caract_curvas).
%  2) Etapa 14 (MFDFA h(+1) + Higuchi) y etapa 15 (LOSO/LOPO, voto por
%     registro, permutaciones de todas las configuraciones, baselines y
%     orden vs nivel) sobre esas ventanas.
%  Los días con registros legados no entran al control de orden vs nivel:
%  su marca "Created on" es la fecha de guardado/conversión, no la de grabación.
%  -> resultados_asr/legados/
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
dirL = fullfile(cfg.rutas.resultados, 'legados');
dirC = fullfile(dirL, 'cache');
if ~isfolder(dirC), mkdir(dirC); end
nw = usar_paralelo(cfg);
registrar(cfg, '==== Corrida con registros legados (14 KMI / 14 FIS) ====');

Cand = candidatos_rango(512, cfg);
baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for b = {'Re64a2', 'Re128a4'}
    c = Cand(strcmp({Cand.nombre}, b{1}));
    baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<SAGROW>
end
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
ses = ses(logical(ses.completa), :);
dirS = fullfile(cfg.rutas.resultados, 'senales');

%% 1) Ventanas contiguas ampliadas
for L = [512 1000]
    pre = sprintf('W%d', L);
    destino = fullfile(dirC, ['w_' pre '.mat']);
    if isfile(destino), continue, end
    t0 = tic;
    S = load(fullfile(cfg.rutas.resultados, 'complejidad', 'cache', ['w_' pre '.mat']));
    W = S.W;
    nuevas = ses(~ismember(ses.sesion, unique(W.sesion)), :);
    [~, lay] = caract_registro(cfg, L, baseH);
    extractor = @(V) single(caract_curvas(V, cfg, lay));
    rech = struct('id', 'R0', 'umbral_z', Inf, 'min_canales', 1);
    partes = cell(height(nuevas), 1); sid = nuevas.sesion; cnd = nuevas.condicion;
    parfor (s = 1:height(nuevas), nw)
        partes{s} = ventanas_tiempo_sesion(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, baseH, L, cfg, ...
            rech, extractor, max(L, 512));
    end
    Wn = vertcat(partes{:});
    Wn.feat = single(matriz_caracteristicas(Wn));
    W = [W; Wn];
    save(destino, 'W', '-v7.3');
    registrar(cfg, '  %s: +%d sesiones (%s), %d ventanas nuevas, %d en total (%.1f min)', pre, height(nuevas), ...
        strjoin(nuevas.sesion, ', '), height(Wn), height(W), toc(t0)/60);
end

%% 2) Etapa 14 y etapa 15 sobre las ventanas ampliadas
cfg.mfh.carpeta = fullfile('legados', 'mfdfa_higuchi');
cfg.mfh.dir_w   = fullfile('legados', 'cache');
etapa14_mfdfa_higuchi(cfg);
cfg.revision.carpeta = fullfile('legados', 'revision');
cfg.revision.dir_w   = fullfile('legados', 'cache');
cfg.revision.perm_todas = true;
cfg.revision.dias_excluidos_orden = {'Px009_M1', 'Px010_M1'};
etapa15_revision(cfg);
registrar(cfg, '==== Corrida con registros legados terminada ====');
