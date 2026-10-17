function etapa9_fractales(cfg)
%ETAPA9_FRACTALES  ¿Qué información añaden los algoritmos fractales de Niveles3?
%   Amplía la matriz de características con los 8 algoritmos fractales de
%   Niveles3 (RS con destendencia, Higuchi, DFA, pendiente del variograma,
%   su promedio, HO, HRS64 y HV; ver fractales_ventana.m) y mide su aporte
%   frente al par calibrado de Hurst (resultados/calibracion/mejor_conjunto.txt).
%
%   A) Matrices de características (ventanas de 512 muestras):
%      1. KMI 3 clases intra-sesión: MISMAS ventanas que los ARFF de la
%         etapa 3 (anclas y balanceo leídos de los *_meta.csv), como en la
%         etapa 7 -> arff/fractales/intra/*.arff
%      2. 10 vs 40 entre sesiones: ventanas contiguas de todo el registro
%         (ventanas_tiempo_sesion) -> fractales/ventanas_fractales.mat
%   B) Evaluación escalonada de conjuntos de familias:
%      base | cada algoritmo solo | base + cada algoritmo | los 8 de
%      Niveles3 | base + los 8 | selección hacia adelante desde la base.
%      Métrica 1: KMI 3c, CV 5-fold, 13 sesiones x 3 clasificadores
%                 ajustados (sin selección de atributos) + CV anidada de la
%                 elección del conjunto.
%      Métrica 2: generalización 10 vs 40 (evaluar_generalizacion: LOSO KMI
%                 y FIS, Wilcoxon, Px010-M1 invertido, excedente sobre el
%                 control de mitades) con el SVM ajustado (ig8 y sin selección).
%      Métrica 3: matriz de confusión 3c y herencia del reposo (base vs mejor).
%      Permutación (1000) de la LOSO KMI del mejor conjunto si supera a la base.
%   -> resultados/fractales/*.csv, figuras/

dirF  = fullfile(cfg.rutas.resultados, 'fractales');
dirC  = fullfile(dirF, 'cache');
dirFg = fullfile(dirF, 'figuras');
dirA  = fullfile(cfg.rutas.resultados, 'arff', 'fractales', 'intra');
for d = {dirF, dirC, dirFg, dirA}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
L = cfg.ajuste.duracion;
nw = usar_paralelo(cfg);
ch = cfg.canales.nombres;

% Par calibrado (base) a partir de su nombre
nomBase = strtrim(fileread(fullfile(cfg.rutas.resultados, 'calibracion', 'mejor_conjunto.txt')));
Cand = candidatos_rango(L, cfg);
base = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for b = strsplit(nomBase, '+')
    c = Cand(strcmp({Cand.nombre}, b{1}));
    base(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<AGROW>
end
famN3 = cfg.fractales.familias;
extractor = @(V) fractales_ventana(V, cfg, base);
registrar(cfg, 'Etapa 9: algoritmos fractales de Niveles3 (%s) frente a la base %s.', strjoin(famN3, ', '), nomBase);
avance(dirF, 'pasos', 0, 6);

%% A1) ARFF KMI 3c intra (mismas ventanas que la etapa 3)
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
    meta = readtable(fullfile(dirArff, strrep(Ii.archivo{s}, '.arff', '_meta.csv')), 'TextType', 'char');
    [X, nombres] = caracteristicas_intra(fullfile(dirS, [Ii.id{s} '.mat']), meta, extractor, L, cfg);
    y = arrayfun(@(c) sprintf('N%02d', c), meta.clase, 'UniformOutput', false);
    exportar_arff(rutas{s}, Ii.dataset{s}, nombres, X, {'N00','N10','N40'}, y, ...
        {sprintf('Etapa 9: base %s + algoritmos fractales de Niveles3', nomBase)});
    copyfile(fullfile(dirArff, strrep(Ii.archivo{s}, '.arff', '_meta.csv')), strrep(rutas{s}, '.arff', '_meta.csv'));
end
registrar(cfg, '  ARFF 3c intra: %d sesiones (%.1f min).', ns, toc(t0)/60);

%% A2) Ventanas contiguas para la generalización 10 vs 40
archV = fullfile(dirF, 'ventanas_fractales.mat');
if isfile(archV)
    S = load(archV); W = S.W;
else
    t0 = tic;
    ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
    ses = ses(logical(ses.completa), :);
    if cfg.modo_prueba, ses = ses(1:min(4, height(ses)), :); end
    rech = struct('id', 'R0', 'umbral_z', Inf, 'min_canales', 1);
    partes = cell(height(ses), 1);
    sid = ses.sesion; cnd = ses.condicion;
    parfor (s = 1:height(ses), nw)
        partes{s} = ventanas_tiempo_sesion(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, base, L, cfg, rech, extractor);
    end
    W = vertcat(partes{:});
    save(archV, 'W', '-v7.3');
    registrar(cfg, '  ventanas contiguas: %d (%.1f min).', height(W), toc(t0)/60);
end
[~, nombresW] = extractor(zeros(numel(ch), L) + randn(numel(ch), L));
avance(dirF, 'pasos', 1, 6);

%% B1) Conjuntos fijos: KMI 3c intra
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste', 'hiperparametros_elegidos.csv'), 'TextType', 'char', 'Delimiter', ',');
H.Properties.VariableNames = {'modelo', 'especificacion'};
H = H(~strcmp(H.modelo, 'RF (ajustado)'), :);
modelos = arrayfun(@(i) modelo_desde_id(H.especificacion{i}, 'ninguna'), 1:height(H));
etiq = strtok(H.modelo).';
n3 = strjoin(famN3, '+');
conj = [{nomBase}, famN3, strcat(nomBase, '+', famN3), {n3, [nomBase '+' n3]}];
tipo = [{'base'}, repmat({'solo'}, 1, numel(famN3)), repmat({'base+1'}, 1, numel(famN3)), {'Niveles3', 'base+Niveles3'}];
t0 = tic;
cache = fullfile(dirC, strcat('conj_', Ii.id, '.mat'));
parfor (s = 1:ns, nw)
    evaluar_conjuntos_cv(rutas{s}, conj, modelos, etiq, cfg, cache{s});
end
registrar(cfg, '  3c: %d conjuntos x %d sesiones x %d clasificadores (%.1f min).', numel(conj), ns, numel(modelos), toc(t0)/60);
avance(dirF, 'pasos', 2, 6);

%% B2) Selección hacia adelante desde la base (3c)
actual = nomBase;
restantes = famN3;
T = cargar_cache_conjuntos(dirC, 'conj_', Ii.id);
Rs = resumir_conjuntos(T);
accActual = Rs.acc_media(strcmp(Rs.conjunto, actual));
pasos = {0, actual, accActual, ''};
while ~isempty(restantes)
    cand = strcat(actual, '+', restantes);
    parfor (s = 1:ns, nw)
        evaluar_conjuntos_cv(rutas{s}, cand, modelos, etiq, cfg, cache{s});
    end
    T = cargar_cache_conjuntos(dirC, 'conj_', Ii.id);
    Rs = resumir_conjuntos(T);
    [~, ia] = ismember(cand, Rs.conjunto);
    [mejorAcc, j] = max(Rs.acc_media(ia));
    if mejorAcc <= accActual, break, end
    actual = cand{j}; accActual = mejorAcc;
    pasos(end+1, :) = {size(pasos, 1), actual, accActual, restantes{j}}; %#ok<AGROW>
    registrar(cfg, '  adelante paso %d: + %s -> %.2f %%', size(pasos, 1) - 1, restantes{j}, accActual);
    restantes(j) = [];
end
writetable(cell2table(pasos, 'VariableNames', {'paso','conjunto','acc_media','familia_anadida'}), ...
    fullfile(dirF, 'seleccion_adelante.csv'));
if ~ismember(actual, conj)
    conj{end+1} = actual; tipo{end+1} = 'adelante';
end
T = cargar_cache_conjuntos(dirC, 'conj_', Ii.id);
writetable(T, fullfile(dirF, 'fractales_3c_por_sesion.csv'));
Rs = resumir_conjuntos(T);
Rs = Rs(ismember(Rs.conjunto, conj), :);
[~, it] = ismember(Rs.conjunto, conj);
Rs.tipo = tipo(it).';
Rs.n_atributos = numel(ch) * (1 + count(Rs.conjunto, '+'));
Rs.aporte_sobre_base = Rs.acc_media - Rs.acc_media(strcmp(Rs.conjunto, nomBase));
registrar(cfg, '  3c base %s = %.2f %%; mejor %s = %.2f %%.', nomBase, ...
    Rs.acc_media(strcmp(Rs.conjunto, nomBase)), Rs.conjunto{1}, Rs.acc_media(1));
avance(dirF, 'pasos', 3, 6);

%% B3) CV anidada de la elección del conjunto + matriz de confusión / herencia del reposo
mejor3c = Rs.conjunto{1};
porClf = T(strcmp(T.conjunto, mejor3c), :);
[g, cl] = findgroups(porClf.clasificador);
[~, j] = max(splitapply(@mean, porClf.acc, g));
mClf = modelos(strcmp(etiq, cl{j}));
archAn = fullfile(dirF, 'fractales_3c_anidada.csv');
if ~isfile(archAn)
    t0 = tic;
    filas = cell(ns, 1);
    parfor (s = 1:ns, nw)
        filas{s} = cv_anidada_conjuntos(rutas{s}, conj, mClf, cfg);
    end
    A = table(Ii.id, repmat(cl(j), ns, 1), cellfun(@(f) f.acc, filas), ...
        cellfun(@(f) strjoin(f.elegidos, ' | '), filas, 'UniformOutput', false), ...
        'VariableNames', {'sesion','clasificador','acc_anidada','elegidos_por_fold'});
    writetable(A, archAn);
    registrar(cfg, '  CV anidada de la elección (%s): %.2f +- %.2f %% (%.1f min).', cl{j}, ...
        mean(A.acc_anidada), std(A.acc_anidada), toc(t0)/60);
end
A = readtable(archAn, 'TextType', 'char');
Rs.acc_anidada_mejor = nan(height(Rs), 1);
Rs.acc_anidada_mejor(1) = mean(A.acc_anidada);
conf = {};
for cj = unique({nomBase, mejor3c}, 'stable')
    [M, her] = confusion_herencia(rutas, cj{1}, mClf, cfg, nw);
    conf(end+1, :) = {cj{1}, cl{j}, M(1,1), M(1,2), M(1,3), M(2,1), M(2,2), M(2,3), M(3,1), M(3,2), M(3,3), her}; %#ok<AGROW>
    registrar(cfg, '  confusión %s (%s): reposo bien = %.1f %%, herencia del reposo = %.1f %%', ...
        cj{1}, cl{j}, M(1,1), her);
end
writetable(cell2table(conf, 'VariableNames', {'conjunto','clasificador','r00_p00','r00_p10','r00_p40', ...
    'r10_p00','r10_p10','r10_p40','r40_p00','r40_p10','r40_p40','herencia_reposo'}), ...
    fullfile(dirF, 'fractales_confusion_3c.csv'));
writetable(Rs, fullfile(dirF, 'fractales_3c_resumen.csv'));
avance(dirF, 'pasos', 4, 6);

%% B4) Generalización 10 vs 40 (LOSO) por conjunto y modelo
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
specSVM = Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')};
mods = {modelo_desde_id(specSVM), modelo_desde_id(specSVM, 'ninguna')};
etiqM = {'SVM ajustado ig8', 'SVM ajustado sin seleccion'};
trabajos = {};
for c = 1:numel(conj), for m = 1:numel(mods), trabajos(end+1, :) = {c, m}; end, end %#ok<AGROW>
nt = size(trabajos, 1);
archG = arrayfun(@(k) fullfile(dirC, sprintf('gen_%03d.mat', k)), (1:nt).', 'UniformOutput', false);
t0 = tic;
pendG = find(~cellfun(@isfile, archG));
parfor (q = 1:numel(pendG), nw)
    k = pendG(q);
    c = trabajos{k, 1}; m = trabajos{k, 2};
    cols = columnas(nombresW, conj{c});
    Tk = W; Tk.feat = matriz_caracteristicas(W); Tk.feat = Tk.feat(:, cols);
    fila = [table(conj(c), tipo(c), etiqM(m), numel(cols), 'VariableNames', {'conjunto','tipo','modelo','n_atributos'}), ...
        evaluar_generalizacion(Tk, nombresW(cols), mods{m}, cfg)];
    guardar_fila(archG{k}, fila);
end
G = cellfun(@(a) getfield(load(a), 'fila'), archG, 'UniformOutput', false); %#ok<GFLD>
G = vertcat(G{:});
G = sortrows(G, 'loso_kmi', 'descend');
writetable(G, fullfile(dirF, 'fractales_generalizacion.csv'));
registrar(cfg, '  generalización: %d conjuntos x %d modelos (%.1f min). Mejor LOSO KMI: %s (%s) = %.1f %%', ...
    numel(conj), numel(mods), toc(t0)/60, G.conjunto{1}, G.modelo{1}, G.loso_kmi(1));
avance(dirF, 'pasos', 5, 6);

%% B5) Permutación del mejor conjunto (LOSO KMI, modelo ig8) si supera a la base
% Solo con el SVM ig8: sin selección, una LOSO con ~160 atributos tarda
% ~10 min y 1000 permutaciones serían ~1 día (prueba de rendimiento).
G8 = sortrows(G(strcmp(G.modelo, etiqM{1}), :), 'loso_kmi', 'descend');
gb = G8(strcmp(G8.conjunto, nomBase), :);
archP = fullfile(dirF, 'permutacion_fractales.csv');
if ~strcmp(G8.conjunto{1}, nomBase) && G8.loso_kmi(1) > gb.loso_kmi && ~isfile(archP)
    m = 1;
    cols = columnas(nombresW, G8.conjunto{1});
    Tk = W(strcmp(W.condicion, 'KMI'), :);
    Tk.feat = matriz_caracteristicas(Tk); Tk.feat = Tk.feat(:, cols);
    obs = mean(loso_niveles(Tk, nombresW(cols), mods{m}, cfg, []));
    sesC = unique(Tk.sesion);
    np = 1000; lote = 50; nula = zeros(np, 1);
    t0 = tic;
    for b = 1:lote:np
        idx = b:min(np, b + lote - 1);
        parte = zeros(numel(idx), 1);
        parfor (jj = 1:numel(idx), nw)
            rs = RandStream('mt19937ar', 'Seed', 500 + idx(jj));
            voltear = sesC(rand(rs, numel(sesC), 1) < 0.5);
            parte(jj) = mean(loso_niveles(Tk, nombresW(cols), mods{m}, cfg, voltear));
        end
        nula(idx) = parte;
        avance(dirF, 'permutacion', idx(end), np);
    end
    p = (1 + nnz(nula >= obs)) / (1 + np);
    writetable(table(G8.conjunto(1), G8.modelo(1), obs, mean(nula), std(nula), percentil95(nula), p, ...
        'VariableNames', {'conjunto','modelo','loso_kmi','nula_media','nula_de','nula_p95','p'}), archP);
    registrar(cfg, '  permutación %s (%s): LOSO KMI = %.1f %% | nula %.1f +- %.1f (p95 %.1f) | p = %.4f (%.1f min)', ...
        G8.conjunto{1}, G8.modelo{1}, obs, mean(nula), std(nula), percentil95(nula), p, toc(t0)/60);
else
    registrar(cfg, '  sin permutación: el mejor conjunto en LOSO es la base o no la supera.');
end
graficar(Rs, G, nomBase, famN3, dirFg);
avance(dirF, 'pasos', 6, 6);
registrar(cfg, 'Etapa 9 lista: %s', dirF);
end

%% ======================================================================
function [X, nombres] = caracteristicas_intra(rutaSenal, meta, extractor, L, cfg)
% Recalcula las ventanas del ARFF de la etapa 3 y extrae todas las familias.
S = load(rutaSenal);
nombresArch = cellfun(@(r) regexprep(r, '^.*[\\/]|\.mat$', ''), S.archivos, 'UniformOutput', false);
X = [];
for j = 1:numel(S.X)
    filas = find(strcmp(meta.archivo, nombresArch{j}));
    if isempty(filas), continue, end
    A = segmentar_por_estimulos(size(S.X{j}, 2), S.stims{j}, cfg);
    for f = filas.'
        k = find(A.bloque == meta.bloque(f) & A.ancla == meta.ancla(f), 1);
        [x, nombres] = extractor(double(extraer_ventanas(S.X{j}, A.centro(k), L)));
        if isempty(X), X = nan(height(meta), numel(x)); end
        X(f, :) = x;
    end
end
end

function cols = columnas(nombres, conjunto)
% Índices de los atributos '<familia>_<canal>' de las familias del conjunto.
fam = strsplit(conjunto, '+');
pref = cellfun(@(n) n(1:find(n == '_', 1, 'last') - 1), nombres, 'UniformOutput', false);
cols = find(ismember(pref, fam));
end

function [M, her] = confusion_herencia(rutas, conjunto, m, cfg, nw)
% Matriz de confusión (% por fila) y % del reposo predicho con el nivel de su archivo.
partes = cell(numel(rutas), 1);
parfor (s = 1:numel(rutas), nw)
    weka_iniciar(cfg);
    D = weka_cargar_arff(rutas{s});
    attrs = {};
    fams = strsplit(conjunto, '+');
    for q = 1:numel(fams), attrs = [attrs, strcat(fams{q}, '_', cfg.canales.nombres)]; end %#ok<AGROW>
    D = weka_seleccionar_atributos(D, attrs);
    [pred, real] = weka_predicciones_cv(D, m, cfg, cfg.weka.semilla);
    meta = readtable(strrep(rutas{s}, '.arff', '_meta.csv'), 'TextType', 'char');
    partes{s} = [real, pred, meta.nivel_archivo];
end
P = vertcat(partes{:});
M = zeros(3);
for r = 1:3
    for c = 1:3, M(r, c) = 100 * mean(P(P(:,1) == r, 2) == c); end
end
rep = P(P(:,1) == 1, :);
her = 100 * mean((rep(:,3) == 10 & rep(:,2) == 2) | (rep(:,3) == 40 & rep(:,2) == 3));
end

function guardar_fila(ruta, fila)
tmp = [ruta '.tmp.mat']; save(tmp, 'fila'); movefile(tmp, ruta, 'f');
end

function avance(dirF, nombre, hecho, total)
fid = fopen(fullfile(dirF, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end

function graficar(Rs, G, nomBase, famN3, dirFg)
f = figure('Visible', 'off', 'Position', [100 100 1200 450]);
solo = cellfun(@(a) Rs.acc_media(strcmp(Rs.conjunto, a)), famN3);
mas = cellfun(@(a) Rs.acc_media(strcmp(Rs.conjunto, [nomBase '+' a])), famN3);
b0 = Rs.acc_media(strcmp(Rs.conjunto, nomBase));
subplot(1, 2, 1);
bar(categorical(famN3, famN3), [solo(:), mas(:)]); hold on; yline(b0, '--k', nomBase); yline(100/3, ':', 'azar');
ylabel('Exactitud KMI 3c (%)'); legend({'solo', ['base + algoritmo']}, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('KMI 3 clases intra-sesión');
subplot(1, 2, 2);
Gs = G(strcmp(G.modelo, 'SVM ajustado sin seleccion'), :);
v1 = cellfun(@(a) val(Gs, a), famN3); v2 = cellfun(@(a) val(Gs, [nomBase '+' a]), famN3);
bar(categorical(famN3, famN3), [v1(:), v2(:)]); hold on; yline(val(Gs, nomBase), '--k', nomBase); yline(50, ':', 'azar');
ylabel('LOSO KMI 10 vs 40 (%)'); title('Generalización a sesión no vista (SVM sin selección)');
exportgraphics(f, fullfile(dirFg, 'fractales_aporte.png'), 'Resolution', 150);
close(f);
end

function v = val(G, conjunto)
v = G.loso_kmi(strcmp(G.conjunto, conjunto));
if isempty(v), v = NaN; end
v = v(1);
end
