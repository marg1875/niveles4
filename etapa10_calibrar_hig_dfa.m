function etapa10_calibrar_hig_dfa(cfg)
%ETAPA10_CALIBRAR_HIG_DFA  Calibración de Higuchi y DFA, solos y combinados.
%   En la etapa 9, Higuchi (y en menor medida DFA) fueron los únicos
%   algoritmos fractales que añadían algo a Hurst, con los parámetros por
%   defecto de Niveles3. Aquí se calibran como se hizo con Hurst (etapa 7):
%
%   Candidatos (candidatos_hig_dfa): Higuchi con rango de k = kmin..kmax
%   (72) y DFA con rango contiguo de escalas y destendencia lineal (DFA1) o
%   cuadrática (DFA2) (210); referencias: HIGn3, DFAn3 (Niveles3) y las
%   familias de Hurst de cfg.calib_hd.base_hurst. Todos salen de las curvas
%   L(k) y F(n) de cada ventana (curvas_ventana), calculadas una sola vez.
%
%   1) Curvas en las ventanas 3c intra (etapa 3) y en las ventanas
%      contiguas de 512 (generalización) -> ARFF con todos los candidatos.
%   2) Prueba individual de cada candidato con dos criterios:
%        3c   KMI 3 clases intra-sesión (13 sesiones x kNN/SVM/MLP ajustados)
%        LOSO KMI 10 vs 40 dejando una sesión fuera (SVM ajustado ig8)
%   3) Pares dentro de cada familia y cruzados HIG+DFA (finalistas por 3c).
%   4) Combinación con Hurst: selección hacia adelante desde la base de
%      cada criterio (3c: Re16a2+Re64a4; LOSO: Re64a2+Re128a4) sobre los
%      finalistas de HIG y DFA.
%   5) Controles de sesgo: CV anidada (3c) de la elección dentro de cada
%      familia y de la combinación; LOSO anidada de la elección entre los
%      finalistas (elegidos sin mirar la LOSO); métricas completas de
%      generalización (evaluar_generalizacion) y permutación (1000) de los
%      mejores conjuntos por LOSO.
%   -> resultados/calib_hig_dfa/*.csv, figuras/

p = cfg.calib_hd;
dirK = fullfile(cfg.rutas.resultados, 'calib_hig_dfa');
dirC = fullfile(dirK, 'cache');
dirFg = fullfile(dirK, 'figuras');
dirA = fullfile(cfg.rutas.resultados, 'arff', 'calib_hig_dfa', 'intra');
for d = {dirK, dirC, dirFg, dirA}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
L = cfg.ajuste.duracion;
nw = usar_paralelo(cfg);
Cand = candidatos_rango(L, cfg);
baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for b = p.base_hurst
    c = Cand(strcmp({Cand.nombre}, b{1}));
    baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<AGROW>
end
[K, lay] = candidatos_hig_dfa(cfg, baseH);
nomHIG = {K(strcmp({K.tipo}, 'HIG')).nombre};
nomDFA = {K(strcmp({K.tipo}, 'DFA')).nombre};
base3c = 'Re16a2+Re64a4';
baseLOSO = 'Re64a2+Re128a4';
extractor = @(V) curvas_ventana(V, cfg, baseH, lay);
registrar(cfg, 'Etapa 10: calibración de Higuchi (%d) y DFA (%d) + referencias.', numel(nomHIG), numel(nomDFA));
avance(dirK, 'pasos', 0, 8);

%% 1a) Curvas y ARFF 3c intra (mismas ventanas que la etapa 3)
dirArff = fullfile(cfg.rutas.resultados, 'arff', sprintf('d%d', L));
I = readtable(fullfile(dirArff, 'indice.csv'), 'TextType', 'char');
Ii = I(strcmp(I.condicion, cfg.ajuste.condicion) & strcmp(I.esquema, cfg.ajuste.esquema) & ...
    strcmp(I.alcance, 'intra'), :);
if cfg.modo_prueba, Ii = Ii(1:2, :); end
ns = height(Ii);
rutas = fullfile(dirA, strcat(Ii.dataset, '.arff'));
dirS = fullfile(cfg.rutas.resultados, 'senales');
t0 = tic;
pend = find(~cellfun(@isfile, rutas));
parfor (q = 1:numel(pend), nw)
    s = pend(q);
    archMeta = fullfile(dirArff, strrep(Ii.archivo{s}, '.arff', '_meta.csv'));
    meta = readtable(archMeta, 'TextType', 'char');
    M = curvas_intra(fullfile(dirS, [Ii.id{s} '.mat']), meta, extractor, L, cfg);
    [X, nombres] = caracteristicas_desde_curvas(M, K, lay, cfg);
    y = arrayfun(@(c) sprintf('N%02d', c), meta.clase, 'UniformOutput', false);
    exportar_arff(rutas{s}, Ii.dataset{s}, nombres, X, {'N00','N10','N40'}, y, ...
        {'Etapa 10: calibración de Higuchi y DFA (todos los candidatos)'});
    copyfile(archMeta, strrep(rutas{s}, '.arff', '_meta.csv'));
end
registrar(cfg, '  ARFF 3c intra con %d candidatos x %d canales (%.1f min).', numel(K), lay.C, toc(t0)/60);

%% 1b) Curvas en las ventanas contiguas (generalización)
archV = fullfile(dirK, 'ventanas_curvas.mat');
if isfile(archV)
    S = load(archV); W = S.W;
else
    t0 = tic;
    ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
    ses = ses(logical(ses.completa), :);
    rech = struct('id', 'R0', 'umbral_z', Inf, 'min_canales', 1);
    partes = cell(height(ses), 1);
    sid = ses.sesion; cnd = ses.condicion;
    parfor (s = 1:height(ses), nw)
        partes{s} = ventanas_tiempo_sesion(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, baseH, L, cfg, rech, extractor);
    end
    W = vertcat(partes{:});
    [Xw, nombresW] = caracteristicas_desde_curvas(matriz_caracteristicas(W), K, lay, cfg);
    W.feat = Xw;
    save(archV, 'W', 'nombresW', '-v7.3');
    registrar(cfg, '  ventanas contiguas: %d (%.1f min).', height(W), toc(t0)/60);
end
S = load(archV, 'nombresW'); nombresW = S.nombresW;
WK = W(strcmp(W.condicion, 'KMI'), :);
avance(dirK, 'pasos', 1, 8);

%% 2a) Prueba individual 3c (todos los candidatos)
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste', 'hiperparametros_elegidos.csv'), 'TextType', 'char', 'Delimiter', ',');
H.Properties.VariableNames = {'modelo', 'especificacion'};
H = H(~strcmp(H.modelo, 'RF (ajustado)'), :);
modelos = arrayfun(@(i) modelo_desde_id(H.especificacion{i}, 'ninguna'), 1:height(H));
etiq = strtok(H.modelo).';
cache = fullfile(dirC, strcat('c3_', Ii.id, '.mat'));
individuales = [{K.nombre}, {base3c, baseLOSO}];
t0 = tic;
parfor (s = 1:ns, nw)
    evaluar_conjuntos_cv(rutas{s}, individuales, modelos, etiq, cfg, cache{s});
end
registrar(cfg, '  3c individual: %d candidatos x %d sesiones x %d clasificadores (%.1f min).', ...
    numel(individuales), ns, numel(modelos), toc(t0)/60);
avance(dirK, 'pasos', 2, 8);

%% 2b) Prueba individual LOSO KMI (SVM ajustado ig8)
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
mLOSO = modelo_desde_id(Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')});
t0 = tic;
Lind = loso_conjuntos(WK, nombresW, individuales, mLOSO, cfg, nw, dirC, 'loso_', dirK);
registrar(cfg, '  LOSO individual: %d conjuntos (%.1f min).', numel(individuales), toc(t0)/60);
avance(dirK, 'pasos', 3, 8);

%% Tabla individual
T = cargar_cache_conjuntos(dirC, 'c3_', Ii.id);
R = resumir_conjuntos(T);
R = R(ismember(R.conjunto, individuales), :);
[~, ik] = ismember(R.conjunto, {K.nombre});
tipos = repmat({'REF'}, height(R), 1); tipos(ik > 0) = {K(ik(ik > 0)).tipo};
R.familia = tipos;
[~, il] = ismember(R.conjunto, Lind.conjunto);
R.loso_kmi = Lind.loso_kmi(il);
R.loso_ses_sobre50 = Lind.ses_sobre50(il);
R.rango_3c = rango(R.acc_media);
R.rango_loso = rango(R.loso_kmi);
writetable(R, fullfile(dirK, 'individual.csv'));
for fam = {'HIG', 'DFA'}
    Rf = sortrows(R(strcmp(R.familia, fam{1}), :), 'acc_media', 'descend');
    Rl = sortrows(R(strcmp(R.familia, fam{1}), :), 'loso_kmi', 'descend');
    registrar(cfg, '  %s mejor 3c: %s = %.2f %% (LOSO %.1f) | mejor LOSO: %s = %.1f %% (3c %.2f)', fam{1}, ...
        Rf.conjunto{1}, Rf.acc_media(1), Rf.loso_kmi(1), Rl.conjunto{1}, Rl.loso_kmi(1), Rl.acc_media(1));
end
for r = {'HIGn3', 'DFAn3', base3c, baseLOSO}
    i = strcmp(R.conjunto, r{1});
    registrar(cfg, '  referencia %-15s: 3c %.2f %% | LOSO %.1f %%', r{1}, R.acc_media(i), R.loso_kmi(i));
end

%% 3) Pares (finalistas por 3c)
topH = top(R, 'HIG', p.top_familia);
topD = top(R, 'DFA', p.top_familia);
pares = [pares_de(topH), pares_de(topD)];
for a = 1:numel(topH), for b = 1:numel(topD), pares{end+1} = [topH{a} '+' topD{b}]; end, end %#ok<AGROW>
t0 = tic;
parfor (s = 1:ns, nw)
    evaluar_conjuntos_cv(rutas{s}, pares, modelos, etiq, cfg, cache{s});
end
registrar(cfg, '  3c pares: %d (%.1f min).', numel(pares), toc(t0)/60);
avance(dirK, 'pasos', 4, 8);

%% 4) Combinación con Hurst: selección hacia adelante
pool = [topH, topD];
[pasos3c, final3c] = adelante_3c(base3c, pool, rutas, cache, modelos, etiq, cfg, nw, Ii.id, dirC);
registrar(cfg, '  adelante 3c: %s', final3c);
[pasosL, finalL] = adelante_loso(baseLOSO, pool, WK, nombresW, mLOSO, cfg, nw, dirC, dirK);
registrar(cfg, '  adelante LOSO: %s', finalL);
writetable([pasos3c; pasosL], fullfile(dirK, 'seleccion_adelante.csv'));
avance(dirK, 'pasos', 5, 8);

%% Resumen 3c de todo lo evaluado
T = cargar_cache_conjuntos(dirC, 'c3_', Ii.id);
writetable(T, fullfile(dirK, 'c3_por_sesion.csv'));
Rt = resumir_conjuntos(T);
Rt.n_atributos = lay.C * (1 + count(Rt.conjunto, '+'));
Rt.aporte_base3c = Rt.acc_media - Rt.acc_media(strcmp(Rt.conjunto, base3c));
writetable(Rt, fullfile(dirK, 'c3_resumen.csv'));

%% 5a) CV anidada 3c: elección dentro de cada familia y de las combinaciones
porClf = T(strcmp(T.conjunto, final3c), :);
[g, cl] = findgroups(porClf.clasificador);
[~, j] = max(splitapply(@mean, porClf.acc, g));
mClf = modelos(strcmp(etiq, cl{j}));
grupos = struct('nombre', {'HIG', 'DFA', 'combinaciones'}, ...
    'cand', {nomHIG, nomDFA, unique([{base3c}, strcat(base3c, '+', pool), pasos3c.conjunto.', {final3c}, pares], 'stable')});
filas = {};
for gI = 1:numel(grupos)
    arch = fullfile(dirC, ['anidada_' grupos(gI).nombre '.mat']);
    if isfile(arch)
        S = load(arch); A = S.A;
    else
        t0 = tic; A = cell(ns, 1); cand = grupos(gI).cand;
        parfor (s = 1:ns, nw)
            A{s} = cv_anidada_conjuntos(rutas{s}, cand, mClf, cfg);
        end
        save(arch, 'A');
        registrar(cfg, '  CV anidada 3c (%s, %d candidatos, %s): %.2f %% (%.1f min).', grupos(gI).nombre, ...
            numel(cand), cl{j}, mean(cellfun(@(a) a.acc, A)), toc(t0)/60);
    end
    for s = 1:ns
        filas(end+1, :) = {grupos(gI).nombre, Ii.id{s}, cl{j}, A{s}.acc, strjoin(A{s}.elegidos, ' | ')}; %#ok<AGROW>
    end
end
writetable(cell2table(filas, 'VariableNames', {'grupo','sesion','clasificador','acc_anidada','elegidos_por_fold'}), ...
    fullfile(dirK, 'c3_anidada.csv'));
avance(dirK, 'pasos', 6, 8);

%% 5b) LOSO anidada sobre finalistas elegidos SIN mirar la LOSO (por 3c)
candL = unique([{baseLOSO, base3c, 'HIGn3', 'DFAn3'}, topH, topD, strcat(baseLOSO, '+', [topH, topD])], 'stable');
archAL = fullfile(dirK, 'loso_anidada.csv');
if ~isfile(archAL)
    t0 = tic;
    [accA, eleg] = loso_anidado(WK, nombresW, candL, mLOSO, cfg, nw);
    writetable(table(accA, eleg, 'VariableNames', {'acc','elegido'}), archAL);
    registrar(cfg, '  LOSO anidada (%d candidatos): %.1f %% (%d/%d sesiones > 50) (%.1f min).', numel(candL), ...
        mean(accA), nnz(accA > 50), numel(accA), toc(t0)/60);
end

%% 5c) Generalización completa de los finalistas y permutación
Lt = sortrows(cargar_loso(dirC, 'loso_'), 'loso_kmi', 'descend');
mejorH_L = Lt.conjunto{find(startsWith(Lt.conjunto, 'HIG') & ~contains(Lt.conjunto, '+'), 1)};
mejorD_L = Lt.conjunto{find(startsWith(Lt.conjunto, 'DFA') & ~contains(Lt.conjunto, '+'), 1)};
finalistas = unique({baseLOSO, base3c, 'HIGn3', 'DFAn3', topH{1}, topD{1}, mejorH_L, mejorD_L, ...
    [topH{1} '+' topD{1}], final3c, finalL}, 'stable');
G = generalizacion(finalistas, W, nombresW, mLOSO, cfg, nw, dirC);
writetable(G, fullfile(dirK, 'generalizacion_finalistas.csv'));
avance(dirK, 'pasos', 7, 8);
perm = {};
for cj = unique({finalL, mejorH_L, mejorD_L}, 'stable')
    archP = fullfile(dirC, ['perm_' strrep(cj{1}, '+', '_') '.mat']);
    if ~isfile(archP)
        t0 = tic;
        [obs, nula] = permutacion(WK, nombresW, cj{1}, mLOSO, cfg, nw, p.n_permutaciones, dirK);
        save(archP, 'obs', 'nula');
        registrar(cfg, '  permutación %s: LOSO = %.1f %% | nula %.1f +- %.1f | p = %.4f (%.1f min)', cj{1}, ...
            obs, mean(nula), std(nula), (1 + nnz(nula >= obs)) / (1 + numel(nula)), toc(t0)/60);
    end
    S = load(archP);
    perm(end+1, :) = {cj{1}, S.obs, mean(S.nula), std(S.nula), percentil95(S.nula), ...
        (1 + nnz(S.nula >= S.obs)) / (1 + numel(S.nula))}; %#ok<AGROW>
end
writetable(cell2table(perm, 'VariableNames', {'conjunto','loso_kmi','nula_media','nula_de','nula_p95','p'}), ...
    fullfile(dirK, 'permutacion.csv'));
graficar(R, dirFg);
avance(dirK, 'pasos', 8, 8);
registrar(cfg, 'Etapa 10 lista: %s', dirK);
end

%% ======================================================================
function M = curvas_intra(rutaSenal, meta, extractor, L, cfg)
% Curvas de las ventanas del ARFF de la etapa 3 (mismas anclas).
S = load(rutaSenal);
nombresArch = cellfun(@(r) regexprep(r, '^.*[\\/]|\.mat$', ''), S.archivos, 'UniformOutput', false);
M = [];
for j = 1:numel(S.X)
    filas = find(strcmp(meta.archivo, nombresArch{j}));
    if isempty(filas), continue, end
    A = segmentar_por_estimulos(size(S.X{j}, 2), S.stims{j}, cfg);
    for f = filas.'
        k = find(A.bloque == meta.bloque(f) & A.ancla == meta.ancla(f), 1);
        x = extractor(double(extraer_ventanas(S.X{j}, A.centro(k), L)));
        if isempty(M), M = nan(height(meta), numel(x)); end
        M(f, :) = x;
    end
end
end

function c = top(R, fam, n)
Rf = sortrows(R(strcmp(R.familia, fam), :), 'acc_media', 'descend');
c = Rf.conjunto(1:min(n, height(Rf))).';
end

function P = pares_de(c)
P = {};
for a = 1:numel(c), for b = a+1:numel(c), P{end+1} = [c{a} '+' c{b}]; end, end %#ok<AGROW>
end

function cols = columnas(nombres, conjunto)
pref = cellfun(@(n) n(1:find(n == '_', 1, 'last') - 1), nombres, 'UniformOutput', false);
cols = find(ismember(pref, strsplit(conjunto, '+')));
end

function Lt = loso_conjuntos(WK, nombresW, conjuntos, m, cfg, nw, dirC, prefijo, dirK)
% LOSO KMI de cada conjunto, con caché por conjunto (reanudable).
arch = fullfile(dirC, strcat(prefijo, strrep(conjuntos, '+', '_'), '.mat'));
pend = find(~cellfun(@isfile, arch));
X = matriz_caracteristicas(WK); WKm = removevars(WK, 'feat');
hechos = numel(conjuntos) - numel(pend);
for b = 1:7*4:numel(pend)                         % lotes para el avance
    idx = pend(b:min(numel(pend), b + 7*4 - 1));
    parfor (q = 1:numel(idx), nw)
        k = idx(q);
        cols = columnas(nombresW, conjuntos{k});
        Tk = WKm; Tk.feat = X(:, cols);
        acc = loso_niveles(Tk, nombresW(cols), m, cfg, []);
        fila = table(conjuntos(k), mean(acc), nnz(acc > 50), numel(acc), ...
            'VariableNames', {'conjunto','loso_kmi','ses_sobre50','n_ses'});
        guardar(arch{k}, fila);
    end
    hechos = hechos + numel(idx);
    avance(dirK, 'loso', hechos, numel(conjuntos));
end
Lt = cargar_loso(dirC, prefijo, arch);
end

function Lt = cargar_loso(dirC, prefijo, arch)
if nargin < 3
    d = dir(fullfile(dirC, [prefijo '*.mat']));
    d = d(~endsWith({d.name}, '.tmp.mat'));
    arch = fullfile(dirC, {d.name});
end
F = cellfun(@(a) getfield(load(a), 'fila'), arch, 'UniformOutput', false); %#ok<GFLD>
Lt = vertcat(F{:});
end

function [pasos, actual] = adelante_3c(base, pool, rutas, cache, modelos, etiq, cfg, nw, ids, dirC)
% Selección hacia adelante por exactitud 3c media.
actual = base; restantes = pool;
Rs = resumir_conjuntos(cargar_cache_conjuntos(dirC, 'c3_', ids));
acc = Rs.acc_media(strcmp(Rs.conjunto, actual));
pasos = table({'3c'}, 0, {actual}, acc, {''}, 'VariableNames', {'criterio','paso','conjunto','acc','anadido'});
while ~isempty(restantes)
    cand = strcat(actual, '+', restantes);
    parfor (s = 1:numel(rutas), nw)
        evaluar_conjuntos_cv(rutas{s}, cand, modelos, etiq, cfg, cache{s});
    end
    Rs = resumir_conjuntos(cargar_cache_conjuntos(dirC, 'c3_', ids));
    [~, ia] = ismember(cand, Rs.conjunto);
    [a, j] = max(Rs.acc_media(ia));
    if a <= acc, break, end
    actual = cand{j}; acc = a;
    pasos = [pasos; table({'3c'}, height(pasos), {actual}, acc, restantes(j), ...
        'VariableNames', pasos.Properties.VariableNames)]; %#ok<AGROW>
    restantes(j) = [];
end
end

function [pasos, actual] = adelante_loso(base, pool, WK, nombresW, m, cfg, nw, dirC, dirK)
% Selección hacia adelante por LOSO KMI (ig8).
actual = base; restantes = pool;
L0 = loso_conjuntos(WK, nombresW, {actual}, m, cfg, nw, dirC, 'loso_', dirK);
acc = L0.loso_kmi;
pasos = table({'LOSO'}, 0, {actual}, acc, {''}, 'VariableNames', {'criterio','paso','conjunto','acc','anadido'});
while ~isempty(restantes)
    cand = strcat(actual, '+', restantes);
    Lc = loso_conjuntos(WK, nombresW, cand, m, cfg, nw, dirC, 'loso_', dirK);
    [a, j] = max(Lc.loso_kmi);
    if a <= acc, break, end
    actual = cand{j}; acc = a;
    pasos = [pasos; table({'LOSO'}, height(pasos), {actual}, acc, restantes(j), ...
        'VariableNames', pasos.Properties.VariableNames)]; %#ok<AGROW>
    restantes(j) = [];
end
end

function r = rango(x)
% Puesto (1 = mayor) sin Statistics Toolbox; empates con el mismo puesto.
[~, ~, u] = unique(-x);
r = u;
end

function G = generalizacion(conj, W, nombresW, m, cfg, nw, dirC)
arch = fullfile(dirC, strcat('gen_', strrep(conj, '+', '_'), '.mat'));
X = matriz_caracteristicas(W); Wm = removevars(W, 'feat');
pend = find(~cellfun(@isfile, arch));
parfor (q = 1:numel(pend), nw)
    k = pend(q);
    cols = columnas(nombresW, conj{k});
    Tk = Wm; Tk.feat = X(:, cols);
    fila = [table(conj(k), numel(cols), 'VariableNames', {'conjunto','n_atributos'}), ...
        evaluar_generalizacion(Tk, nombresW(cols), m, cfg)];
    guardar(arch{k}, fila);
end
F = cellfun(@(a) getfield(load(a), 'fila'), arch, 'UniformOutput', false); %#ok<GFLD>
G = sortrows(vertcat(F{:}), 'loso_kmi', 'descend');
end

function [obs, nula] = permutacion(WK, nombresW, conjunto, m, cfg, nw, np, dirK)
cols = columnas(nombresW, conjunto);
Tk = removevars(WK, 'feat'); X = matriz_caracteristicas(WK); Tk.feat = X(:, cols);
obs = mean(loso_niveles(Tk, nombresW(cols), m, cfg, []));
sesC = unique(Tk.sesion);
nula = zeros(np, 1);
lote = 50;
for b = 1:lote:np
    idx = b:min(np, b + lote - 1);
    parte = zeros(numel(idx), 1);
    parfor (jj = 1:numel(idx), nw)
        rs = RandStream('mt19937ar', 'Seed', 500 + idx(jj));
        voltear = sesC(rand(rs, numel(sesC), 1) < 0.5);
        parte(jj) = mean(loso_niveles(Tk, nombresW(cols), m, cfg, voltear));
    end
    nula(idx) = parte;
    avance(dirK, ['perm_' strrep(conjunto, '+', '_')], idx(end), np);
end
end

function guardar(ruta, fila)
tmp = [ruta '.tmp.mat']; save(tmp, 'fila'); movefile(tmp, ruta, 'f');
end

function avance(dirK, nombre, hecho, total)
fid = fopen(fullfile(dirK, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end

function graficar(R, dirFg)
% Mapas de calor: exactitud 3c y LOSO por rango de parámetros.
f = figure('Visible', 'off', 'Position', [100 100 1500 900]);
Rh = R(strcmp(R.familia, 'HIG'), :);
tok = regexp(Rh.conjunto, 'HIG(\d+)a(\d+)', 'tokens', 'once');
a = cellfun(@(t) str2double(t{1}), tok); b = cellfun(@(t) str2double(t{2}), tok);
[ua, ~, ia] = unique(a); [ub, ~, ib] = unique(b);
crit = {'acc_media', 'KMI 3c intra (%)'; 'loso_kmi', 'LOSO KMI 10 vs 40 (%)'};
for q = 1:2
    Mq = nan(numel(ua), numel(ub));
    Mq(sub2ind(size(Mq), ia, ib)) = Rh.(crit{q,1});
    subplot(2, 3, q);
    h = heatmap(string(ub), string(ua), round(Mq, 1)); h.XLabel = 'k_{max}'; h.YLabel = 'k_{min}';
    h.Title = ['Higuchi - ' crit{q,2}];
end
for o = 1:2
    letra = 'lc';
    Rd = R(startsWith(R.conjunto, ['DFA' letra(o)]), :);
    tok = regexp(Rd.conjunto, 'DFA[lc](\d+)a(\d+)', 'tokens', 'once');
    a = cellfun(@(t) str2double(t{1}), tok); b = cellfun(@(t) str2double(t{2}), tok);
    [ua, ~, ia] = unique(a); [ub, ~, ib] = unique(b);
    for q = 1:2
        Mq = nan(numel(ua), numel(ub));
        Mq(sub2ind(size(Mq), ia, ib)) = Rd.(crit{q,1});
        subplot(2, 3, 2 + (o-1)*2 + q);
        h = heatmap(string(ub), string(ua), round(Mq, 1)); h.XLabel = 'n_{max}'; h.YLabel = 'n_{min}';
        h.Title = sprintf('DFA%d - %s', o, crit{q,2});
    end
end
exportgraphics(f, fullfile(dirFg, 'calibracion_hig_dfa.png'), 'Resolution', 130);
close(f);
end
