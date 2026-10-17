function etapa11_niveles3_ganadores(cfg)
%ETAPA11_NIVELES3_GANADORES  Conjuntos de caracterizadores ganadores de Niveles3.
%   Niveles3 (INFORME_FRACTAL.md, §5.2-5.3) obtuvo su mejor resultado
%   fractal con "Básicos per-channel" (16 canales x {RS, Higuchi, DFA,
%   Variograma} = 64 atributos) y comparó "Todos 7 métodos" (112) y
%   "Martinez" (HO + HRS64 + HV, 48), con LogisticRegression, SVM RBF
%   (C = 12, gamma = 0.02) y RandomForest (300 árboles), CV intra-sujeto.
%   Aquí se reproducen esos conjuntos y clasificadores en Weka, sobre los
%   mismos datos y controles de este pipeline, junto con las referencias:
%     Re16a2+Re64a4    Hurst calibrado (3c)       Re64a2+Re128a4   Hurst (LOSO)
%     HIG3a8+DFAc4a24  Higuchi + DFA calibrados (etapa 10)
%     RS+HIG3a8+DFAc4a24+VAR  Básicos con Higuchi y DFA calibrados
%   Evaluación:
%     a) KMI 3c intra-sesión y 2c (10 vs 40, sin reposo) intra-sesión,
%        CV 5-fold, 13 sesiones, por mes; 6 clasificadores (3 ajustados de
%        este pipeline + 3 de Niveles3). Los datasets ya están balanceados
%        (balanceo del reposo de la etapa 3), por eso no se usa SMOTE.
%     b) Matriz de confusión y herencia del reposo del mejor 3c.
%     c) Generalización 10 vs 40 (evaluar_generalizacion) con SVM ig8,
%        LogReg y SVM RBF de Niveles3.
%     d) Permutación secuencial (Besag-Clifford) de los 2 mejores LOSO.
%   Usa los ARFF/ventanas de las etapas 9 y 10 (no recalcula características).
%   -> resultados/niveles3_ganadores/*.csv

dirG = fullfile(cfg.rutas.resultados, 'niveles3_ganadores');
dirC = fullfile(dirG, 'cache');
dirA3 = fullfile(cfg.rutas.resultados, 'arff', 'niveles3_ganadores', 'intra3c');
dirA2 = fullfile(cfg.rutas.resultados, 'arff', 'niveles3_ganadores', 'intra2c');
for d = {dirG, dirC, dirA3, dirA2}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
nw = usar_paralelo(cfg);
ch = cfg.canales.nombres;
conj = {'RS+HIG+DFA+VAR', 'RS+HIG+DFA+VAR+HO+HRS64+HV', 'HO+HRS64+HV', ...
        'RS+HIG3a8+DFAc4a24+VAR', 'HIG3a8+DFAc4a24', 'Re16a2+Re64a4', 'Re64a2+Re128a4'};
etiqConj = {'N3 Básicos', 'N3 Todos 7', 'N3 Martinez', 'Básicos calibrados', ...
            'HIG+DFA calibrados', 'Hurst 3c', 'Hurst LOSO'};
fam = unique(strsplit(strjoin(conj, '+'), '+'), 'stable');
registrar(cfg, 'Etapa 11: conjuntos ganadores de Niveles3 (%d conjuntos, %d workers).', numel(conj), nw);
avance(dirG, 'pasos', 0, 5);

%% 1) ARFF combinados (etapas 9 y 10) 3c y 2c
I9 = dir(fullfile(cfg.rutas.resultados, 'arff', 'fractales', 'intra', '*.arff'));
ids = erase({I9.name}, '.arff');
r3 = fullfile(dirA3, strcat(ids, '.arff')).';
r2 = fullfile(dirA2, strcat(ids, '.arff')).';
pend = find(~cellfun(@isfile, r3) | ~cellfun(@isfile, r2));
dir9 = fullfile(cfg.rutas.resultados, 'arff', 'fractales', 'intra');
dir10 = fullfile(cfg.rutas.resultados, 'arff', 'calib_hig_dfa', 'intra');
parfor (q = 1:numel(pend), nw)
    s = pend(q);
    weka_iniciar(cfg);
    [X9, y9, n9] = weka_a_matriz(weka_cargar_arff(fullfile(dir9, [ids{s} '.arff'])));
    [X10, y10, n10] = weka_a_matriz(weka_cargar_arff(fullfile(dir10, [ids{s} '.arff'])));
    assert(isequal(y9, y10), 'etapa11: las ventanas de las etapas 9 y 10 no coinciden (%s)', ids{s});
    [X, nombres] = unir(X9, n9, X10, n10, fam);
    exportar_arff(r3{s}, ids{s}, nombres, X, {'N00','N10','N40'}, y9, {'Etapa 11 (3c)'});
    t = ~strcmp(y9, 'N00');
    exportar_arff(r2{s}, ids{s}, nombres, X(t, :), {'N10','N40'}, y9(t), {'Etapa 11 (2c: sin reposo)'});
    copyfile(fullfile(dir9, [ids{s} '_meta.csv']), strrep(r3{s}, '.arff', '_meta.csv'));
end
avance(dirG, 'pasos', 1, 5);

%% 2) Clasificadores: ajustados de este pipeline + ganadores de Niveles3
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste', 'hiperparametros_elegidos.csv'), 'TextType', 'char', 'Delimiter', ',');
H.Properties.VariableNames = {'modelo', 'especificacion'};
H = H(~strcmp(H.modelo, 'RF (ajustado)'), :);
modelos = arrayfun(@(i) modelo_desde_id(H.especificacion{i}, 'ninguna'), 1:height(H));
etiq = strtok(H.modelo).';
modelos = [modelos, especificar_modelo('LOG', struct('ridge', 0.5)), ...
    especificar_modelo('SVM', struct('kernel', 'rbf', 'C', 12, 'gamma', 0.02, 'filtro', 1)), ...
    especificar_modelo('RF', struct('arboles', 300))];
etiq = [etiq, {'N3-LogReg', 'N3-SVMrbf', 'N3-RF300'}];
t0 = tic;
for esq = {'3c', '2c'}
    rutas = r3; if strcmp(esq{1}, '2c'), rutas = r2; end
    cache = fullfile(dirC, strcat(['c' esq{1} '_'], ids, '.mat'));
    parfor (s = 1:numel(rutas), nw)
        evaluar_conjuntos_cv(rutas{s}, conj, modelos, etiq, cfg, cache{s});
    end
    T = cargar_cache_conjuntos(dirC, ['c' esq{1} '_'], ids);
    T.esquema = repmat(esq, height(T), 1);
    T.mes = cellfun(@(s) str2double(regexp(s, '_M(\d)', 'tokens', 'once')), T.sesion);
    T.etiqueta = etiqConj(cellfun(@(c) find(strcmp(conj, c)), T.conjunto)).';
    writetable(T, fullfile(dirG, ['intra_' esq{1} '_por_sesion.csv']));
    tabla_intra(T, conj, etiqConj, etiq, fullfile(dirG, ['intra_' esq{1} '_resumen.csv']));
end
registrar(cfg, '  intra 3c y 2c: %d conjuntos x %d clasificadores x %d sesiones (%.1f min).', ...
    numel(conj), numel(modelos), numel(ids), toc(t0)/60);
avance(dirG, 'pasos', 2, 5);

%% 3) Matriz de confusión 3c y herencia del reposo (mejor 3c y referencias)
T3 = readtable(fullfile(dirG, 'intra_3c_por_sesion.csv'), 'TextType', 'char');
[g, cj, cl] = findgroups(T3.conjunto, T3.clasificador);
media = splitapply(@mean, T3.acc, g);
[~, j] = max(media);
filas = {};
for c = unique([cj(j), {'RS+HIG+DFA+VAR', 'Re16a2+Re64a4'}], 'stable')
    k = find(strcmp(etiq, cl{j}));
    [M, her] = confusion_herencia(r3, c{1}, modelos(k), cfg, nw);
    filas(end+1, :) = {c{1}, cl{j}, M(1,1), M(2,2), M(3,3), mean(diag(M)), her}; %#ok<AGROW>
    registrar(cfg, '  confusión %s (%s): reposo %.1f | 10 %.1f | 40 %.1f | herencia %.1f %%', ...
        c{1}, cl{j}, M(1,1), M(2,2), M(3,3), her);
end
writetable(cell2table(filas, 'VariableNames', {'conjunto','clasificador','reposo','n10','n40','media_diag','herencia_reposo'}), ...
    fullfile(dirG, 'confusion_3c.csv'));
avance(dirG, 'pasos', 3, 5);

%% 4) Generalización 10 vs 40 (ventanas contiguas de las etapas 9 y 10)
S9 = load(fullfile(cfg.rutas.resultados, 'fractales', 'ventanas_fractales.mat'));
S10 = load(fullfile(cfg.rutas.resultados, 'calib_hig_dfa', 'ventanas_curvas.mat'));
base9 = struct('nombre', {'Re16a2','Re64a4'}, 'modo', {'ecuacion13','ecuacion13'}, 'm_max', {16, 64}, 'm_min', {2, 4});
[~, n9] = fractales_ventana(randn(numel(ch), cfg.ajuste.duracion), cfg, base9);
W = S9.W;
assert(isequal(W.sesion, S10.W.sesion) && max(abs(W.t_min - S10.W.t_min)) < 1e-9, ...
    'etapa11: las ventanas contiguas de las etapas 9 y 10 no coinciden');
[Xw, nombresW] = unir(matriz_caracteristicas(W), n9, matriz_caracteristicas(S10.W), S10.nombresW, fam);
W.feat = Xw;
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
modsG = [modelo_desde_id(Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')}), modelos(end-2:end-1)];
etiqG = {'SVM ajustado ig8', 'N3-LogReg', 'N3-SVMrbf'};
trab = {};
for c = 1:numel(conj), for m = 1:numel(modsG), trab(end+1, :) = {c, m}; end, end %#ok<AGROW>
arch = arrayfun(@(k) fullfile(dirC, sprintf('gen_%02d.mat', k)), (1:size(trab, 1)).', 'UniformOutput', false);
pend = find(~cellfun(@isfile, arch));
Wm = removevars(W, 'feat');
t0 = tic;
parfor (q = 1:numel(pend), nw)
    k = pend(q); c = trab{k, 1}; m = trab{k, 2};
    cols = columnas(nombresW, conj{c});
    Tk = Wm; Tk.feat = Xw(:, cols);
    fila = [table(conj(c), etiqConj(c), etiqG(m), numel(cols), ...
        'VariableNames', {'conjunto','etiqueta','modelo','n_atributos'}), ...
        evaluar_generalizacion(Tk, nombresW(cols), modsG(m), cfg)];
    guardar(arch{k}, fila);
end
F = cellfun(@(a) getfield(load(a), 'fila'), arch, 'UniformOutput', false); %#ok<GFLD>
G = sortrows(vertcat(F{:}), 'loso_kmi', 'descend');
writetable(G, fullfile(dirG, 'generalizacion.csv'));
registrar(cfg, '  generalización: %d trabajos (%.1f min). Mejor: %s / %s = %.1f %%', size(trab, 1), ...
    toc(t0)/60, G.etiqueta{1}, G.modelo{1}, G.loso_kmi(1));
avance(dirG, 'pasos', 4, 5);

%% 5) Permutación secuencial de los 2 mejores por LOSO KMI
WK = W(strcmp(W.condicion, 'KMI'), :);
P = {};
for i = 1:min(2, height(G))
    archP = fullfile(dirC, sprintf('perm_%d.mat', i));
    m = modsG(strcmp(etiqG, G.modelo{i}));
    if ~isfile(archP)
        t0 = tic;
        cols = columnas(nombresW, G.conjunto{i});
        Tk = removevars(WK, 'feat'); Tk.feat = matriz_caracteristicas(WK); Tk.feat = Tk.feat(:, cols);
        [obs, nula, p] = permutacion_secuencial(Tk, nombresW(cols), m, cfg, nw, 1000, 20, ...
            fullfile(dirG, sprintf('perm_%d.avance.txt', i)));
        save(archP, 'obs', 'nula', 'p');
        registrar(cfg, '  permutación %s / %s: LOSO %.1f %% | nula %.1f +- %.1f (%d perm.) | p = %.4f (%.1f min)', ...
            G.etiqueta{i}, G.modelo{i}, obs, mean(nula), std(nula), numel(nula), p, toc(t0)/60);
    end
    S = load(archP);
    P(end+1, :) = {G.conjunto{i}, G.etiqueta{i}, G.modelo{i}, S.obs, mean(S.nula), std(S.nula), numel(S.nula), S.p}; %#ok<AGROW>
end
writetable(cell2table(P, 'VariableNames', {'conjunto','etiqueta','modelo','loso_kmi','nula_media','nula_de','n_perm','p'}), ...
    fullfile(dirG, 'permutacion.csv'));
avance(dirG, 'pasos', 5, 5);
registrar(cfg, 'Etapa 11 lista: %s', dirG);
end

%% ======================================================================
function [X, nombres] = unir(Xa, na, Xb, nb, fam)
% Columnas de las familias 'fam' tomadas de A (etapa 9) o, si no, de B (etapa 10).
pa = prefijos(na); pb = prefijos(nb);
X = []; nombres = {};
for f = 1:numel(fam)
    ia = find(strcmp(pa, fam{f}));
    if ~isempty(ia)
        X = [X, Xa(:, ia)]; nombres = [nombres, na(ia)]; %#ok<AGROW>
    else
        ib = find(strcmp(pb, fam{f}));
        assert(~isempty(ib), 'etapa11: familia %s no encontrada', fam{f});
        X = [X, Xb(:, ib)]; nombres = [nombres, nb(ib)]; %#ok<AGROW>
    end
end
end

function p = prefijos(n)
p = cellfun(@(x) x(1:find(x == '_', 1, 'last') - 1), n, 'UniformOutput', false);
end

function cols = columnas(nombres, conjunto)
cols = find(ismember(prefijos(nombres), strsplit(conjunto, '+')));
end

function tabla_intra(T, conj, etiqConj, etiq, ruta)
% Exactitud media por conjunto x clasificador (todas las sesiones y por mes).
filas = {};
for c = 1:numel(conj)
    for k = 1:numel(etiq)
        s = strcmp(T.conjunto, conj{c}) & strcmp(T.clasificador, etiq{k});
        fila = {etiqConj{c}, conj{c}, etiq{k}, mean(T.acc(s)), std(T.acc(s)), mean(T.kappa(s))};
        for m = [1 3 6], fila{end+1} = mean(T.acc(s & T.mes == m)); end %#ok<AGROW>
        filas(end+1, :) = fila; %#ok<AGROW>
    end
end
R = cell2table(filas, 'VariableNames', {'etiqueta','conjunto','clasificador','acc','acc_de','kappa','acc_mes1','acc_mes3','acc_mes6'});
writetable(sortrows(R, 'acc', 'descend'), ruta);
end

function [M, her] = confusion_herencia(rutas, conjunto, m, cfg, nw)
partes = cell(numel(rutas), 1);
parfor (s = 1:numel(rutas), nw)
    weka_iniciar(cfg);
    D = weka_cargar_arff(rutas{s});
    fams = strsplit(conjunto, '+'); attrs = {};
    for q = 1:numel(fams), attrs = [attrs, strcat(fams{q}, '_', cfg.canales.nombres)]; end %#ok<AGROW>
    [pred, real] = weka_predicciones_cv(weka_seleccionar_atributos(D, attrs), m, cfg, cfg.weka.semilla);
    meta = readtable(strrep(rutas{s}, '.arff', '_meta.csv'), 'TextType', 'char');
    partes{s} = [real, pred, meta.nivel_archivo];
end
P = vertcat(partes{:});
M = zeros(3);
for r = 1:3, for c = 1:3, M(r, c) = 100 * mean(P(P(:,1) == r, 2) == c); end, end
rep = P(P(:,1) == 1, :);
her = 100 * mean((rep(:,3) == 10 & rep(:,2) == 2) | (rep(:,3) == 40 & rep(:,2) == 3));
end

function guardar(ruta, fila)
tmp = [ruta '.tmp.mat']; save(tmp, 'fila'); movefile(tmp, ruta, 'f');
end

function avance(dirG, nombre, hecho, total)
fid = fopen(fullfile(dirG, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end
