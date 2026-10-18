function etapa13_complejidad(cfg)
%ETAPA13_COMPLEJIDAD  Búsqueda amplia de caracterizadores: todas las familias, todos los parámetros.
%   29 familias (caract_registro): escalamiento (Higuchi, DFA1-3 +/- traslape,
%   DMA, MFDFA h(q) y ancho, varianza agregada, momentos absolutos, Hurst
%   por wavelets, exponente 1/f, box-counting), dimensiones directas (Katz,
%   Petrosian, Sevcik, Hjorth, line length) y entropías/complejidad
%   (SampEn, FuzzyEn, PE, WPE, Lempel-Ziv, MSE), cada una con su rejilla
%   de parámetros, en ventanas de cfg.complejidad.L muestras.
%
%   Fase 1  criterio 3c: KMI 3 clases intra-sesión (13 sesiones x kNN/SVM/MLP
%           ajustados) para TODOS los candidatos (anclas de la etapa 3; para
%           L > 512 solo las anclas 1 y 3, sin traslape).
%   Fase 2  LOSO KMI 10 vs 40 (SVM ig8) para los cfg.complejidad.top_loso
%           mejores de cada familia y L por el criterio 3c (+ referencias).
%   Fase 3  rankings por familia; anidado estricto (elegir por 3c sin la
%           sesión de prueba, medir su LOSO en ella) global, por familia y
%           por L; selección hacia adelante ENTRE familias, global y anidada
%           (toda la selección se repite sin cada sesión de prueba).
%   Fase 4  generalización completa (evaluar_generalizacion) y permutación
%           secuencial de los ganadores.
%   Todo es reanudable (caché por L, sesión y bloque de candidatos).
%   -> resultados/complejidad/*.csv

p = cfg.complejidad;
dirX = fullfile(cfg.rutas.resultados, p.carpeta);
dirC = fullfile(dirX, 'cache');
for d = {dirX, dirC}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
nw = usar_paralelo(cfg);
Cand = candidatos_rango(512, cfg);
baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for b = {'Re64a2', 'Re128a4'}
    c = Cand(strcmp({Cand.nombre}, b{1}));
    baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<AGROW>
end
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste', 'hiperparametros_elegidos.csv'), 'TextType', 'char', 'Delimiter', ',');
H.Properties.VariableNames = {'modelo', 'especificacion'};
H = H(~strcmp(H.modelo, 'RF (ajustado)'), :);
mods3c = arrayfun(@(i) modelo_desde_id(H.especificacion{i}, 'ninguna'), 1:height(H));
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
mSVM = modelo_desde_id(Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')});
dirArff = fullfile(cfg.rutas.resultados, 'arff', 'd512');
I = readtable(fullfile(dirArff, 'indice.csv'), 'TextType', 'char');
Ii = I(strcmp(I.condicion, 'KMI') & strcmp(I.esquema, '3c') & strcmp(I.alcance, 'intra'), :);
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
ses = ses(logical(ses.completa), :);
dirS = fullfile(cfg.rutas.resultados, 'senales');
nL = numel(p.L);
registrar(cfg, '==== Etapa 13: complejidad, %d familias, L = %s, %d workers ====', 29, mat2str(p.L), nw);
avance(dirX, 'pasos', 0, 2 * nL + 3);

%% Fases 1 y 2 por longitud de ventana
R = table();
for iL = 1:nL
    L = p.L(iL); pre = sprintf('W%d', L);
    [K, lay] = caract_registro(cfg, L, baseH);
    extractor = @(V) single(caract_curvas(V, cfg, lay));
    % curvas 3c y contiguas
    arch3 = fullfile(dirC, ['c3_' pre '.mat']);
    if ~isfile(arch3)
        t0 = tic; D3 = cell(height(Ii), 1);
        parfor (s = 1:height(Ii), nw)
            meta = readtable(fullfile(dirArff, strrep(Ii.archivo{s}, '.arff', '_meta.csv')), 'TextType', 'char');
            if L > 512, meta = meta(ismember(meta.ancla, [1 3]), :); end
            D3{s} = struct('id', Ii.id{s}, 'sesion', mapear(Ii.id{s}), ...
                'M', curvas_anclas(fullfile(dirS, [Ii.id{s} '.mat']), meta, extractor, L, cfg), ...
                'y', {arrayfun(@(c) sprintf('N%02d', c), meta.clase, 'UniformOutput', false)});
        end
        save(arch3, 'D3', '-v7.3');
        registrar(cfg, '  %s: curvas 3c (%d ventanas, %d candidatos) en %.1f min', pre, ...
            sum(cellfun(@(d) numel(d.y), D3)), numel(K), toc(t0)/60);
    end
    archW = fullfile(dirC, ['w_' pre '.mat']);
    if ~isfile(archW)
        t0 = tic; partes = cell(height(ses), 1); sid = ses.sesion; cnd = ses.condicion;
        rech = struct('id', 'R0', 'umbral_z', Inf, 'min_canales', 1);
        parfor (s = 1:height(ses), nw)
            partes{s} = ventanas_tiempo_sesion(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, baseH, L, cfg, ...
                rech, extractor, max(L, 512));
        end
        W = vertcat(partes{:});
        W.feat = single(matriz_caracteristicas(W));
        save(archW, 'W', '-v7.3');
        registrar(cfg, '  %s: %d ventanas contiguas en %.1f min', pre, height(W), toc(t0)/60);
    end
    % fase 1: criterio 3c (caché por sesión)
    S = load(arch3); D3 = S.D3;
    archA = cellfun(@(d) fullfile(dirC, sprintf('A_%s_%s.mat', pre, d.id)), D3, 'UniformOutput', false);
    pend = find(~cellfun(@isfile, archA));
    if ~isempty(pend)
        t0 = tic;
        parfor (q = 1:numel(pend), nw)
            s = pend(q);
            A = criterio_3c(D3{s}, K, lay, cfg, mods3c);
            guardar(archA{s}, A);
        end
        registrar(cfg, '  %s: fase 1 (3c) %d candidatos x %d sesiones (%.1f min)', pre, numel(K), numel(D3), toc(t0)/60);
    end
    A = vertcat(cellfun(@(a) getfield(load(a), 'R'), archA, 'UniformOutput', false)); %#ok<GFLD>
    A = vertcat(A{:});
    avance(dirX, 'pasos', 2 * iL - 1, 2 * nL + 3);
    % fase 2: LOSO de los mejores por familia
    [g, cj] = findgroups(A.conjunto);
    r3 = table(cj, splitapply(@mean, A.acc3c, g), 'VariableNames', {'conjunto','acc3c'});
    [~, ik] = ismember(r3.conjunto, {K.nombre});
    r3.familia = {K(ik).familia}.';
    sel = {};
    for fam = unique(r3.familia).'
        rf = sortrows(r3(strcmp(r3.familia, fam{1}), :), 'acc3c', 'descend');
        sel = [sel; rf.conjunto(1:min(p.top_loso, height(rf)))]; %#ok<AGROW>
    end
    B = loso_candidatos(sel, K, lay, archW, cfg, mSVM, nw, dirC, dirX, pre);
    R = [R; resumen(A, B, K, L)]; %#ok<AGROW>
    avance(dirX, 'pasos', 2 * iL, 2 * nL + 3);
end
writetable(R, fullfile(dirX, 'ranking_candidatos.csv'));

%% Fase 3a: rankings por familia y anidado estricto
Rf = ranking_familias(R);
writetable(Rf, fullfile(dirX, 'ranking_familias.csv'));
for i = 1:height(Rf)
    registrar(cfg, '  %2d. %-8s mejor 3c %-28s %.1f %% (LOSO %.1f) | mejor LOSO %-28s %.1f %% (3c %.1f)', i, ...
        Rf.familia{i}, Rf.mejor_3c{i}, Rf.acc3c(i), Rf.loso_de_mejor_3c(i), Rf.mejor_loso{i}, Rf.loso(i), Rf.acc3c_de_mejor_loso(i));
end
[Aall, Ball] = cargar_todo(dirC, p.L, Ii);
pools = [{'todos'}, unique(R.familia).', arrayfun(@(L) sprintf('W%d', L), p.L, 'UniformOutput', false)];
N = anidado(Aall, Ball, R, pools, p, baseH, cfg, mSVM, nw, dirC, dirX);
writetable(N, fullfile(dirX, 'anidado.csv'));
for i = 1:height(N)
    if strcmp(N.pool{i}, 'todos') || startsWith(N.pool{i}, 'W') || N.acc(i) >= 58
        registrar(cfg, '  anidado %-8s: LOSO %.1f %% (%d/%d > 50, p_W = %.3f) | más elegido %s', N.pool{i}, ...
            N.acc(i), N.sobre50(i), N.n(i), N.p_wilcoxon(i), N.mas_elegido{i});
    end
end
avance(dirX, 'pasos', 2 * nL + 1, 2 * nL + 3);

%% Fase 3b: selección hacia adelante entre familias (global y anidada), por L
F = table();
for iL = 1:nL
    F = [F; adelante(p.L(iL), R, p, baseH, cfg, mods3c, mSVM, nw, dirC, dirX)]; %#ok<AGROW>
end
writetable(F, fullfile(dirX, 'adelante.csv'));
avance(dirX, 'pasos', 2 * nL + 2, 2 * nL + 3);

%% Fase 4: validación final
fin = finalistas(R, Rf, F);
fin = fin(1:min(numel(fin), p.max_finalistas));
V = validacion(fin, p, baseH, cfg, mSVM, nw, dirC, dirX);
writetable(V, fullfile(dirX, 'validacion_final.csv'));
avance(dirX, 'pasos', 2 * nL + 3, 2 * nL + 3);
registrar(cfg, '==== Etapa 13 lista: %s ====', dirX);
end

%% ======================================================================
function M = curvas_anclas(rutaSenal, meta, extractor, L, cfg)
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
        if isempty(M), M = nan(height(meta), numel(x), 'single'); end
        M(f, :) = x;
    end
end
end

function s = mapear(id3c)
t = regexp(id3c, '(Px\d+)_M(\d)', 'tokens', 'once');
s = sprintf('%s_M%s_KMI', t{1}, t{2});
end

function R = criterio_3c(D, K, lay, cfg, mods)
% Exactitud CV 5-fold (3 clasificadores) de cada candidato en una sesión 3c.
weka_iniciar(cfg);
X = caract_desde_curvas(double(D.M), K, lay, cfg);
C = lay.C;
nom = arrayfun(@(i) sprintf('a%d', i), 1:C, 'UniformOutput', false);
acc = nan(numel(K), numel(mods));
for k = 1:numel(K)
    Xk = X(:, (k-1)*C + (1:C));
    if mean(isnan(Xk), 'all') > 0.5, continue, end
    Dk = weka_instancias(D.id, nom, Xk, {'N00','N10','N40'}, D.y);
    for m = 1:numel(mods)
        acc(k, m) = weka_validar(Dk, weka_clasificador_param(mods(m), cfg), cfg).acc;
    end
end
R = table(repmat({D.sesion}, numel(K), 1), {K.nombre}.', mean(acc, 2), acc, ...
    'VariableNames', {'sesion','conjunto','acc3c','acc3c_clf'});
end

function B = loso_candidatos(sel, K, lay, archW, cfg, m, nw, dirC, dirX, pre)
% LOSO por sesión de los candidatos 'sel' (caché por candidato).
arch = fullfile(dirC, strcat('L_', sel, '.mat'));
pend = find(~cellfun(@isfile, arch));
if ~isempty(pend)
    t0 = tic;
    [~, ik] = ismember(sel, {K.nombre});
    lote = 6 * max(nw, 1);
    for b0 = 1:lote:numel(pend)
        grupo = pend(b0:min(numel(pend), b0 + lote - 1));
        parfor (q = 1:numel(grupo), nw)
            i = grupo(q);
            [Mk, meta] = cache_ventanas_kmi(archW);       % cada worker lo lee una vez
            T = meta; T.feat = caract_desde_curvas(Mk, K(ik(i)), lay, cfg);
            acc = loso_por_sesion(T, strcat(sel{i}, '_', cfg.canales.nombres), m, cfg);
            fila = table(repmat(sel(i), numel(acc), 1), unique(meta.sesion, 'stable'), acc, ...
                'VariableNames', {'conjunto','sesion','loso'});
            guardar(arch{i}, fila);
        end
        avance(dirX, ['loso_' pre], b0 + numel(grupo) - 1, numel(pend));
    end
    registrar(cfg, '  %s: fase 2 (LOSO) %d candidatos (%.1f min)', pre, numel(pend), toc(t0)/60);
end
B = cellfun(@(a) getfield(load(a), 'R'), arch, 'UniformOutput', false); %#ok<GFLD>
B = vertcat(B{:});
end

function Rs = resumen(A, B, K, L)
[g, cj] = findgroups(A.conjunto);
Rs = table(cj, repmat(L, numel(cj), 1), splitapply(@(x) mean(x, 'omitnan'), A.acc3c, g), ...
    'VariableNames', {'conjunto','L','acc3c'});
clf = splitapply(@(x) mean(x, 1, 'omitnan'), A.acc3c_clf, g);
Rs.acc3c_kNN = clf(:, 1); Rs.acc3c_SVM = clf(:, 2); Rs.acc3c_MLP = clf(:, 3);
[~, ik] = ismember(Rs.conjunto, {K.nombre});
Rs.familia = {K(ik).familia}.';
[g2, cj2] = findgroups(B.conjunto);
lo = splitapply(@(x) mean(x, 'omitnan'), B.loso, g2);
s50 = splitapply(@(x) nnz(x > 50), B.loso, g2);
[~, ib] = ismember(Rs.conjunto, cj2);
Rs.loso = nan(height(Rs), 1); Rs.sobre50 = Rs.loso;
Rs.loso(ib > 0) = lo(ib(ib > 0)); Rs.sobre50(ib > 0) = s50(ib(ib > 0));
end

function Rf = ranking_familias(R)
filas = {};
for fam = unique(R.familia).'
    r = R(strcmp(R.familia, fam{1}), :);
    [a3, i3] = max(r.acc3c);
    rl = r(isfinite(r.loso), :);
    [lo, il] = max(rl.loso);
    filas(end+1, :) = {fam{1}, height(r), r.conjunto{i3}, a3, r.loso(i3), rl.conjunto{il}, lo, rl.acc3c(il)}; %#ok<AGROW>
end
Rf = cell2table(filas, 'VariableNames', {'familia','n_candidatos','mejor_3c','acc3c','loso_de_mejor_3c', ...
    'mejor_loso','loso','acc3c_de_mejor_loso'});
Rf = sortrows(Rf, 'acc3c', 'descend');
end

function [A, B] = cargar_todo(dirC, Ls, Ii)
A = table();
for L = Ls
    for s = 1:height(Ii)
        S = load(fullfile(dirC, sprintf('A_W%d_%s.mat', L, Ii.id{s}))); A = [A; S.R]; %#ok<AGROW>
    end
end
d = dir(fullfile(dirC, 'L_*.mat')); d = d(~endsWith({d.name}, '.tmp.mat'));
B = cellfun(@(a) getfield(load(fullfile(dirC, a)), 'R'), {d.name}, 'UniformOutput', false); %#ok<GFLD>
B = vertcat(B{:});
end

function T = anidado(A, B, R, pools, p, baseH, cfg, m, nw, dirC, dirX)
% Para cada sesión s: candidato con mayor 3c medio sin s (dentro del pool)
% y su LOSO en s. Si ese candidato no tiene LOSO, se calcula.
sesL = unique(B.sesion, 'stable');
falt = {};
eleg = struct();
for ip = 1:numel(pools)
    cand = pool_de(R, pools{ip});
    Ac = A(ismember(A.conjunto, cand), :);
    e = cell(numel(sesL), 1);
    for s = 1:numel(sesL)
        T0 = Ac(~strcmp(Ac.sesion, sesL{s}), :);
        [g, cj] = findgroups(T0.conjunto);
        [~, j] = max(splitapply(@(x) mean(x, 'omitnan'), T0.acc3c, g));
        e{s} = cj{j};
    end
    eleg.(sprintf('p%d', ip)) = e;
    falt = [falt; e(~ismember(e, B.conjunto))]; %#ok<AGROW>
end
falt = unique(falt);
if ~isempty(falt)
    for L = unique(cellfun(@(c) str2double(regexp(c, '^W(\d+)', 'tokens', 'once')), falt)).'
        [K, lay] = caract_registro(cfg, L, baseH);
        fl = falt(startsWith(falt, sprintf('W%d', L)) & ~startsWith(falt, sprintf('W%d0', L)));
        B = [B; loso_candidatos(fl, K, lay, fullfile(dirC, sprintf('w_W%d.mat', L)), cfg, m, nw, dirC, dirX, sprintf('W%d', L))]; %#ok<AGROW>
    end
end
filas = {};
for ip = 1:numel(pools)
    e = eleg.(sprintf('p%d', ip));
    acc = nan(numel(sesL), 1);
    for s = 1:numel(sesL)
        r = B.loso(strcmp(B.conjunto, e{s}) & strcmp(B.sesion, sesL{s}));
        if ~isempty(r), acc(s) = r(1); end
    end
    ok = isfinite(acc);
    if nnz(ok) < 3, continue, end
    Wt = wilcoxon_holm([acc(ok), 50 * ones(nnz(ok), 1)], {'a', 'b'});
    [u, ~, ic] = unique(e(ok)); [~, im] = max(accumarray(ic, 1));
    filas(end+1, :) = {pools{ip}, mean(acc(ok)), nnz(acc(ok) > 50), nnz(ok), Wt.p, u{im}, strjoin(e(ok).', ' | ')}; %#ok<AGROW>
end
T = cell2table(filas, 'VariableNames', {'pool','acc','sobre50','n','p_wilcoxon','mas_elegido','elegidos'});
T = sortrows(T, 'acc', 'descend');
end

function c = pool_de(R, pool)
if strcmp(pool, 'todos'), c = R.conjunto;
elseif startsWith(pool, 'W') && all(isstrprop(pool(2:end), 'digit')), c = R.conjunto(R.L == str2double(pool(2:end)));
else, c = R.conjunto(strcmp(R.familia, pool));
end
end

function F = adelante(L, R, p, baseH, cfg, mods, m, nw, dirC, dirX)
% Selección hacia adelante entre familias por 3c: global (13 sesiones) y
% anidada (repetida sin cada sesión de prueba; su LOSO en esa sesión).
arch = fullfile(dirC, sprintf('adelante_W%d.mat', L));
if isfile(arch), S = load(arch); F = S.F; return, end
t0 = tic;
pre = sprintf('W%d', L);
[K, lay] = caract_registro(cfg, L, baseH);
S = load(fullfile(dirC, ['c3_' pre '.mat'])); D3 = S.D3;
Rl = R(R.L == L, :);
% pool: el mejor candidato de cada familia por 3c (global, solo para precalcular columnas)
fams = unique(Rl.familia);
mejores = cellfun(@(f) mejor_de(Rl, f), fams, 'UniformOutput', false);
sesL = cellfun(@(d) d.sesion, D3, 'UniformOutput', false);
A3 = cargar_A(dirC, pre, D3);
% columnas candidatas: los 3 mejores de cada familia (para que el pool anidado pueda variar)
top3 = {};
for f = fams.'
    rf = sortrows(Rl(strcmp(Rl.familia, f{1}), :), 'acc3c', 'descend');
    top3 = [top3; rf.conjunto(1:min(3, height(rf)))]; %#ok<AGROW>
end
[~, ik] = ismember(top3, {K.nombre});
Xs = cellfun(@(d) caract_desde_curvas(double(d.M), K(ik), lay, cfg), D3, 'UniformOutput', false);
nombres = top3;
% global (q = 0) y anidado (una selección completa sin cada sesión de prueba), en paralelo
nS = numel(sesL);
todo = cell(nS + 1, 1);
parfor (q = 1:nS + 1, nw)
    entr = 1:nS; if q > 1, entr = setdiff(entr, q - 1); end
    todo{q} = forward(entr, A3, top3, Xs, D3, mods, cfg, p);
end
sel.global = todo{1};
res = todo(2:end);
% LOSO de los conjuntos elegidos
SW = load(fullfile(dirC, ['w_' pre '.mat']));
WK = SW.W(strcmp(SW.W.condicion, 'KMI'), :);
[~, ikW] = ismember(top3, {K.nombre});
XW = caract_desde_curvas(double(WK.feat), K(ikW), lay, cfg);
meta = removevars(WK, 'feat');
sesW = unique(meta.sesion, 'stable');
C = lay.C;
colsDe = @(conj) cell2mat(arrayfun(@(i) (i-1)*C + (1:C), find(ismember(nombres, strsplit(conj, '+'))).', 'UniformOutput', false));
T = meta; T.feat = XW(:, colsDe(sel.global));
lg = loso_por_sesion(T, nombres_cols(nombres, sel.global, cfg), m, cfg);
la = nan(numel(sesL), 1);
parfor (s = 1:numel(sesL), nw)
    Ts = meta; Ts.feat = XW(:, colsDe(res{s}));
    a = loso_por_sesion(Ts, nombres_cols(nombres, res{s}, cfg), m, cfg, sesL(s));
    la(s) = a(strcmp(sesW, sesL{s}));
end
ok = isfinite(la);
Wt = wilcoxon_holm([la(ok), 50 * ones(nnz(ok), 1)], {'a', 'b'});
[u, ~, ic] = unique(res(ok)); [~, im] = max(accumarray(ic, 1));
F = table(L, {sel.global}, mean(A3.acc3c(ismember(A3.conjunto, strsplit(sel.global, '+')))), ...
    mean(lg, 'omitnan'), nnz(lg > 50), mean(la(ok)), nnz(la(ok) > 50), nnz(ok), Wt.p, u(im), {strjoin(res(ok).', ' | ')}, ...
    'VariableNames', {'L','global','acc3c_familias','loso_global','sobre50_global','loso_anidado', ...
    'sobre50_anidado','n','p_wilcoxon_anidado','mas_elegido_anidado','elegidos_anidado'});
save(arch, 'F');
registrar(cfg, '  adelante %s: global %s (LOSO %.1f %%) | anidado LOSO %.1f %% (%d/%d > 50, p_W = %.3f) (%.1f min)', ...
    pre, sel.global, F.loso_global, F.loso_anidado, F.sobre50_anidado, F.n, F.p_wilcoxon_anidado, toc(t0)/60);
avance(dirX, ['adelante_' pre], 1, 1);
end

function n = nombres_cols(nombres, conj, cfg)
c = nombres(ismember(nombres, strsplit(conj, '+')));
n = {};
for i = 1:numel(c), n = [n, strcat(c{i}, '_', cfg.canales.nombres)]; end %#ok<AGROW>
end

function A3 = cargar_A(dirC, pre, D3)
A3 = table();
for s = 1:numel(D3)
    S = load(fullfile(dirC, sprintf('A_%s_%s.mat', pre, D3{s}.id))); A3 = [A3; S.R]; %#ok<AGROW>
end
end

function c = mejor_de(R, fam)
r = R(strcmp(R.familia, fam), :);
[~, i] = max(r.acc3c); c = r.conjunto{i};
end

function conj = forward(entr, A3, top3, Xs, D3, mods, cfg, p)
% Hacia adelante con sesiones 'entr': pool = mejor candidato de cada familia
% (3c medio en 'entr'), se agregan familias mientras suba el 3c medio
% (hasta cfg.complejidad.adelante_max familias) entre las 12 mejores.
sesE = cellfun(@(d) d.sesion, D3(entr), 'UniformOutput', false);
At = A3(ismember(A3.sesion, sesE) & ismember(A3.conjunto, top3), :);
[g, cj] = findgroups(At.conjunto);
m3 = splitapply(@mean, At.acc3c, g);
fam = regexprep(cj, '^W\d+', '');
fam = cellfun(@(f) familia_de(f), fam, 'UniformOutput', false);
[uf, ~, gf] = unique(fam);
pool = cell(numel(uf), 1); pm = zeros(numel(uf), 1);
for f = 1:numel(uf)
    i = find(gf == f); [pm(f), j] = max(m3(i)); pool{f} = cj{i(j)};
end
[~, o] = sort(pm, 'descend');
pool = pool(o(1:min(p.adelante_pool, numel(o))));
actual = {pool{1}}; best = pm(o(1));
restantes = pool(2:end);
C = numel(cfg.canales.nombres);
while numel(actual) < p.adelante_max && ~isempty(restantes)
    v = zeros(numel(restantes), 1);
    for r = 1:numel(restantes)
        conj = [actual, restantes(r)];
        cols = cell2mat(cellfun(@(c) (find(strcmp(top3, c)) - 1) * C + (1:C), conj, 'UniformOutput', false));
        a = zeros(numel(entr), 1);
        for q = 1:numel(entr)
            s = entr(q);
            D = weka_instancias('x', arrayfun(@(i) sprintf('a%d', i), 1:numel(cols), 'UniformOutput', false), ...
                Xs{s}(:, cols), {'N00','N10','N40'}, D3{s}.y);
            a(q) = mean(arrayfun(@(mm) weka_validar(D, weka_clasificador_param(mods(mm), cfg), cfg).acc, 1:numel(mods)));
        end
        v(r) = mean(a);
    end
    [vm, j] = max(v);
    if vm <= best, break, end
    best = vm; actual = [actual, restantes(j)]; restantes(j) = []; %#ok<AGROW>
end
conj = strjoin(actual, '+');
end

function f = familia_de(nombre)
% Familia a partir del nombre sin prefijo W<L> (mismo criterio que caract_registro).
pref = {'HIGn3','REF'; 'DFAn3','REF'; 'Re','REF'; 'HIG','HIG'; 'DFAl','DFA1'; 'DFAc','DFA2'; 'DFAt','DFA3'; ...
    'DMA','DMA'; 'MFw','MFancho'; 'MFhm1','MFhm1'; 'MFhm3','MFhm3'; 'MFhm5','MFhm5'; 'MFhp1','MFhp1'; ...
    'MFhp3','MFhp3'; 'MFhp5','MFhp5'; 'AGV','AGV'; 'ABM','ABM'; 'WV','WAV'; 'PSD','PSD'; 'BOX','BOX'; ...
    'KATZ','KATZ'; 'PET','PET'; 'SEV','SEV'; 'HJ','HJORTH'; 'LL','LL'; 'SE','SAMPEN'; 'FE','FUZZYEN'; ...
    'WPE','WPE'; 'PE','PE'; 'LZ','LZC'; 'MSE','MSE'};
f = 'OTRO';
for i = 1:size(pref, 1)
    if startsWith(nombre, pref{i, 1}), f = pref{i, 2}; return, end
end
end

function fin = finalistas(R, Rf, F)
% Mejor por 3c y por LOSO de las 8 primeras familias, conjuntos hacia
% adelante y referencias.
fin = {};
for i = 1:min(8, height(Rf)), fin = [fin, Rf.mejor_3c(i), Rf.mejor_loso(i)]; end %#ok<AGROW>
fin = [fin, F.global.', strcat('W512', {'Re64a2+'}, 'W512Re128a4')];
r = sortrows(R(isfinite(R.loso), :), 'loso', 'descend');
fin = [fin, r.conjunto(1:3).'];
fin = unique(fin, 'stable');
end

function V = validacion(fin, p, baseH, cfg, m, nw, dirC, dirX)
% Métricas completas (KMI y FIS) y permutación secuencial de cada finalista.
filas = cell(numel(fin), 1);
for i = 1:numel(fin)
    archF = fullfile(dirC, sprintf('final_%02d.mat', i));
    if isfile(archF), S = load(archF); filas{i} = S.fila; continue, end
    t0 = tic;
    partes = strsplit(fin{i}, '+');
    L = str2double(regexp(partes{1}, '^W(\d+)', 'tokens', 'once'));
    if any(~startsWith(partes, sprintf('W%d', L))), continue, end   % (todas las del conjunto son de la misma L)
    [K, lay] = caract_registro(cfg, L, baseH);
    [~, ik] = ismember(partes, {K.nombre});
    S = load(fullfile(dirC, sprintf('w_W%d.mat', L)));
    X = caract_desde_curvas(double(S.W.feat), K(ik), lay, cfg);
    T = removevars(S.W, 'feat'); T.feat = X;
    nom = {}; for c = partes, nom = [nom, strcat(c{1}, '_', cfg.canales.nombres)]; end %#ok<AGROW>
    G = evaluar_generalizacion(T, nom, m, cfg);
    TK = T(strcmp(T.condicion, 'KMI'), :);
    [obs, nula, pp] = permutacion_secuencial(TK, nom, m, cfg, nw, p.n_perm, 20, ...
        fullfile(dirX, sprintf('perm_%02d.avance.txt', i)));
    fila = [table(fin(i), numel(nom), obs, mean(nula), numel(nula), pp, ...
        'VariableNames', {'conjunto','n_atributos','loso_perm','nula_media','n_perm','p_perm'}), G];
    save(archF, 'fila'); filas{i} = fila;
    registrar(cfg, '  final %-40s LOSO %.1f %% (p = %.4f, %d perm.) | FIS %.1f %% | Px010 %.0f %% (%.1f min)', ...
        fin{i}, obs, pp, numel(nula), G.loso_fis, G.px010_invertido, toc(t0)/60);
end
V = vertcat(filas{:});
if ~isempty(V), V = sortrows(V, 'loso_perm', 'descend'); end
end

function guardar(ruta, R) %#ok<INUSD>
tmp = [ruta '.tmp.mat']; save(tmp, 'R'); movefile(tmp, ruta, 'f');
end

function avance(dirX, nombre, hecho, total)
fid = fopen(fullfile(dirX, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end
