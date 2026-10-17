function mejor = etapa7_calibrar_caracterizadores(cfg)
%ETAPA7_CALIBRAR_CARACTERIZADORES  Calibración de los parámetros de Hurst (R/S).
%   HO y HRS(p) son la misma operación con distinto rango de escalas: la
%   pendiente de log10<R/S>_m vs log10 m sobre un rango contiguo de tamaños
%   de bloque diádicos [m_min, m_max] (con N = 512: HO = 512..8, HRS64 =
%   512..4). El parámetro a calibrar es, por tanto, el RANGO DE ESCALAS y la
%   forma de promediar R/S (Ec. 13 o rangors original). HV se omite.
%
%   1) Rehace las MISMAS ventanas de los ARFF intra KMI-3c (anclas y
%      balanceo de reposo idénticos, leídos de los *_meta.csv de la etapa 3)
%      sobre las señales limpias y calcula los 56 caracterizadores
%      candidatos (28 rangos x 2 promedios) x 16 canales
%      -> resultados/arff/calibracion/{intra,global}/*.arff
%   2) Evalúa cada candidato (16 atributos) con los 3 clasificadores
%      ajustados en la etapa 6 (kNN, SVM, MLP; sin selección de atributos),
%      CV 5-fold semilla 1, en las 13 sesiones.
%   3) Evalúa los pares de los cfg.calibracion.top_pares mejores (estilo
%      HRS+HO).
%   4) Elige el mejor conjunto: mayor precisión media (sesiones x 3
%      clasificadores); empate -> menos atributos.
%   5) Control de sesgo: CV anidada de la ELECCIÓN del caracterizador (para
%      el clasificador con mejor resultado): en cada fold externo se elige el
%      candidato con CV interna usando solo el entrenamiento.
%   -> resultados/calibracion/*.csv, figuras y mejor_conjunto.txt
%   Si cfg.calibracion.ajustar_despues, ejecuta la etapa 6 con el conjunto
%   calibrado (salida en resultados/ajuste_calibrado).

sufijo = '';
if cfg.modo_prueba, sufijo = '_prueba'; cfg.calibracion.ajustar_despues = false; end
dirK = fullfile(cfg.rutas.resultados, ['calibracion' sufijo]);
dirKC = fullfile(dirK, 'cache');
dirKF = fullfile(dirK, 'figuras');
dirArffK = fullfile(cfg.rutas.resultados, 'arff', ['calibracion' sufijo]);
for d = {dirK, dirKC, dirKF, fullfile(dirArffK, 'intra'), fullfile(dirArffK, 'global')}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
L = cfg.ajuste.duracion;
C = candidatos_rango(L, cfg);
nw = usar_paralelo(cfg);
registrar(cfg, 'Etapa 7: calibración de caracterizadores (%d candidatos x %d canales).', ...
    numel(C), numel(cfg.canales.nombres));

%% 1) ARFF de calibración
dirArff = fullfile(cfg.rutas.resultados, 'arff', sprintf('d%d', L));
I = readtable(fullfile(dirArff, 'indice.csv'), 'TextType', 'char');
Ii = I(strcmp(I.condicion, cfg.ajuste.condicion) & strcmp(I.esquema, cfg.ajuste.esquema) & ...
    strcmp(I.alcance, 'intra'), :);
if cfg.modo_prueba, Ii = Ii(1:2, :); end
ns = height(Ii);
nombres = {};
for k = 1:numel(C), nombres = [nombres, strcat(C(k).nombre, '_', cfg.canales.nombres)]; end %#ok<AGROW>
if ~isfile(fullfile(dirArffK, 'indice.csv'))
    t0 = tic;
    dirS = fullfile(cfg.rutas.resultados, 'senales');
    partes = cell(ns, 1);
    parfor (s = 1:ns, nw)
        meta = readtable(fullfile(dirArff, strrep(Ii.archivo{s}, '.arff', '_meta.csv')), 'TextType', 'char');
        partes{s} = caracteristicas_calibracion(fullfile(dirS, [Ii.id{s} '.mat']), meta, C, L, cfg);
    end
    filasI = cell(ns + 1, 8);
    clases = {'N00','N10','N40'};
    Xg = []; yg = {}; Mg = table();
    for s = 1:ns
        P = partes{s};
        nom = Ii.dataset{s};
        exportar_arff(fullfile(dirArffK, 'intra', [nom '.arff']), nom, nombres, P.X, clases, P.y, ...
            {'Calibracion de caracterizadores (rangos de escala R/S); ventanas de la etapa 3'});
        writetable(P.meta, fullfile(dirArffK, 'intra', [nom '_meta.csv']));
        filasI(s,:) = {nom, cfg.ajuste.condicion, cfg.ajuste.esquema, 'intra', Ii.id{s}, numel(P.y), ...
            mat2str(cellfun(@(c) nnz(strcmp(P.y, c)), clases)), fullfile('intra', [nom '.arff'])};
        Xg = [Xg; P.X]; yg = [yg; P.y]; Mg = [Mg; P.meta]; %#ok<AGROW>
    end
    nomG = sprintf('%s_%s_GLOBAL', cfg.ajuste.condicion, cfg.ajuste.esquema);
    exportar_arff(fullfile(dirArffK, 'global', [nomG '.arff']), nomG, nombres, Xg, clases, yg, {});
    writetable(Mg, fullfile(dirArffK, 'global', [nomG '_meta.csv']));
    filasI(ns+1,:) = {nomG, cfg.ajuste.condicion, cfg.ajuste.esquema, 'global', 'GLOBAL', numel(yg), ...
        mat2str(cellfun(@(c) nnz(strcmp(yg, c)), clases)), fullfile('global', [nomG '.arff'])};
    writetable(cell2table(filasI, 'VariableNames', {'dataset','condicion','esquema','alcance','id', ...
        'n_instancias','n_por_clase','archivo'}), fullfile(dirArffK, 'indice.csv'));
    registrar(cfg, '  ARFF de calibración listos (%.1f min).', toc(t0)/60);
end
rutas = fullfile(dirArffK, 'intra', strcat(Ii.dataset, '.arff'));

%% 2) Candidatos individuales con los clasificadores ajustados
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste', 'hiperparametros_elegidos.csv'), 'TextType', 'char', ...
    'Delimiter', ',');
H.Properties.VariableNames = {'modelo', 'especificacion'};    % el archivo puede traer BOM
H = H(~strcmp(H.modelo, 'RF (ajustado)'), :);
modelos = arrayfun(@(i) modelo_desde_id(H.especificacion{i}, 'ninguna'), 1:height(H));
etiq = strtok(H.modelo).';
conjuntos = {C.nombre};
t0 = tic;
parfor (s = 1:ns, nw)
    evaluar_conjuntos_cv(rutas{s}, conjuntos, modelos, etiq, cfg, ...
        fullfile(dirKC, ['individual_' Ii.id{s} '.mat']));
end
T1 = cargar_cache_conjuntos(dirKC, 'individual_', Ii.id);
registrar(cfg, '  candidatos individuales evaluados (%.1f min).', toc(t0)/60);

%% 3) Pares de los mejores
R1 = resumir_conjuntos(T1);
top = R1.conjunto(1:min(cfg.calibracion.top_pares, height(R1)));
pares = {};
for a = 1:numel(top)
    for b = a+1:numel(top), pares{end+1} = [top{a} '+' top{b}]; end %#ok<AGROW>
end
t0 = tic;
parfor (s = 1:ns, nw)
    evaluar_conjuntos_cv(rutas{s}, pares, modelos, etiq, cfg, fullfile(dirKC, ['pares_' Ii.id{s} '.mat']));
end
T2 = cargar_cache_conjuntos(dirKC, 'pares_', Ii.id);
registrar(cfg, '  pares evaluados (%.1f min).', toc(t0)/60);

%% 4) Selección
T = [T1; T2];
writetable(T, fullfile(dirK, 'calibracion_por_sesion.csv'));
Rs = resumir_conjuntos(T);
Rs.n_atributos = numel(cfg.canales.nombres) * (1 + count(Rs.conjunto, '+'));
Rs = sortrows(Rs, {'acc_media','n_atributos'}, {'descend','ascend'});
Rs.rango = (1:height(Rs)).';
writetable(Rs, fullfile(dirK, 'calibracion_resumen.csv'));
mejor = Rs.conjunto{1};
ref = {'Re512a8', 'Re512a4', 'Re512a4+Re512a8'};       % HO, HRS64, HRS64+HO del artículo
for r = ref
    i = find(strcmp(Rs.conjunto, r{1}), 1);
    if ~isempty(i)
        registrar(cfg, '  referencia %-16s: %.2f %% (puesto %d de %d)', r{1}, Rs.acc_media(i), i, height(Rs));
    end
end
registrar(cfg, '  MEJOR: %s = %.2f %% (media de %d sesiones x %d clasificadores)', mejor, ...
    Rs.acc_media(1), ns, numel(modelos));
fid = fopen(fullfile(dirK, 'mejor_conjunto.txt'), 'w'); fprintf(fid, '%s\n', mejor); fclose(fid);
graficar_calibracion(Rs, C, dirKF);

%% 5) Control de sesgo: CV anidada de la elección del caracterizador
archAn = fullfile(dirK, 'calibracion_anidada.csv');
if ~isfile(archAn)
    porClf = T(strcmp(T.conjunto, mejor), :);
    [g, cl] = findgroups(porClf.clasificador);
    [~, j] = max(splitapply(@mean, porClf.acc, g));
    m = modelos(strcmp(etiq, cl{j}));
    cand = [conjuntos, pares];
    t0 = tic;
    filas = cell(ns, 1);
    parfor (s = 1:ns, nw)
        filas{s} = cv_anidada_conjuntos(rutas{s}, cand, m, cfg);
    end
    A = table(Ii.id, repmat(cl(j), ns, 1), cellfun(@(f) f.acc, filas), ...
        cellfun(@(f) strjoin(f.elegidos, ' | '), filas, 'UniformOutput', false), ...
        'VariableNames', {'sesion','clasificador','acc_anidada','elegidos_por_fold'});
    writetable(A, archAn);
    registrar(cfg, '  CV anidada de la calibración (%s): %.2f +- %.2f %% (%.1f min).', cl{j}, ...
        mean(A.acc_anidada), std(A.acc_anidada), toc(t0)/60);
end

%% 6) Ajuste de clasificadores con el conjunto calibrado
if cfg.calibracion.ajustar_despues
    c6 = cfg;
    c6.ajuste.dir_arff = dirArffK;
    c6.ajuste.conjunto = mejor;
    c6.ajuste.etiqueta = '_calibrado';
    c6.ajuste.diag_conjuntos = unique([strsplit(mejor, '+'), {mejor}], 'stable');
    etapa6_ajuste_clasificadores(c6);
end
end

%% ======================================================================
function P = caracteristicas_calibracion(rutaSenal, meta, C, L, cfg)
% Recalcula las ventanas del ARFF original y los 56 caracterizadores.
S = load(rutaSenal);
nombresArch = cellfun(@(r) regexprep(r, '^.*[\\/]|\.mat$', ''), S.archivos, 'UniformOutput', false);
ms = floor(L ./ 2.^(0:floor(log2(L))));
ms = ms(ms >= cfg.hurst.bloque_min_abs);
nC = numel(cfg.canales.nombres);
X = zeros(height(meta), numel(C) * nC);
for j = 1:numel(S.X)
    filas = find(strcmp(meta.archivo, nombresArch{j}));
    if isempty(filas), continue, end
    A = segmentar_por_estimulos(size(S.X{j}, 2), S.stims{j}, cfg);
    for f = filas.'
        k = find(A.bloque == meta.bloque(f) & A.ancla == meta.ancla(f), 1);
        V = double(extraer_ventanas(S.X{j}, A.centro(k), L));
        RS = {rs_diadico(V, ms, 'ecuacion13'), rs_diadico(V, ms, 'rangors_original')};
        fila = zeros(1, numel(C) * nC);
        for c = 1:numel(C)
            q = 1 + strcmp(C(c).modo, 'rangors_original');
            idx = C(c).ihi:C(c).ilo;
            fila((c-1)*nC + (1:nC)) = pendiente_hurst(ms(idx), RS{q}(:, idx), cfg.hurst.min_escalas).';
        end
        X(f, :) = fila;
    end
end
P.X = X;
P.y = arrayfun(@(c) sprintf('N%02d', c), meta.clase, 'UniformOutput', false);
P.meta = meta;
end

function graficar_calibracion(Rs, C, dirKF)
f = figure('Visible', 'off', 'Position', [100 100 1100 420]);
modos = {'e', 'Ec. 13 (media de R/S)'; 'o', 'rangors original'};
for q = 1:2
    Cq = C(endsWith(arrayfun(@(c) c.nombre(2), C, 'UniformOutput', false), modos{q,1}));
    mmax = unique([Cq.m_max]); mmin = unique([Cq.m_min]);
    Mq = nan(numel(mmax), numel(mmin));
    for c = Cq
        i = find(strcmp(Rs.conjunto, c.nombre), 1);
        if ~isempty(i), Mq(mmax == c.m_max, mmin == c.m_min) = Rs.acc_media(i); end
    end
    subplot(1, 2, q);
    h = heatmap(string(mmin), string(mmax), round(Mq, 1));
    h.XLabel = 'bloque mínimo m_{min}'; h.YLabel = 'bloque máximo m_{max}';
    h.Title = sprintf('Precisión media (%%) - %s', modos{q,2});
end
exportgraphics(f, fullfile(dirKF, 'calibracion_rangos.png'), 'Resolution', 150);
close(f);
end
