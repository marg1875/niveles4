function etapa6_ajuste_clasificadores(cfg)
%ETAPA6_AJUSTE_CLASIFICADORES  Ajuste de clasificadores para KMI en 3 clases.
%   Datos: ARFF intra-sujeto de la etapa 3 (cfg.ajuste: KMI, 3c, 512
%   muestras, configuración E9, conjunto HRS64+HO+HV) -> 13 sesiones.
%
%   1) Diagnóstico de convergencia: MLP (entrenamiento vs CV por épocas y
%      tasa de aprendizaje, con/sin parada temprana), RF (n.º de árboles) y
%      comparación de familias de Hurst.
%   2) Búsqueda en rejilla sobre la MISMA CV 5-fold (semilla 1), decidida por
%      el usuario: kNN, SVM (poly y RBF), MLP y RF x selección de
%      características dentro de cada fold. Criterio COMÚN a todas las
%      sesiones: mayor precisión media (empate -> mayor kappa medio).
%   3) Ensamble (Vote, promedio de probabilidades) de los 4 mejores.
%   4) Estimación final: CV 5-fold repetida 10 veces para línea base
%      (parámetros del artículo), modelos ajustados y ensamble.
%   5) Control de sesgo: CV anidada (externa 5 x interna 5) para los
%      cfg.ajuste.n_anidada mejores clasificadores ajustados.
%   6) Estadística: binomial por sesión vs azar (1/3), Wilcoxon vs azar y vs
%      línea base (Holm), permutación de etiquetas en el dataset global.
%   Reanudable (cachés en resultados/ajuste/cache). Salidas en resultados/ajuste.

dirA = fullfile(cfg.rutas.resultados, ['ajuste' cfg.ajuste.etiqueta]);
if cfg.modo_prueba, dirA = fullfile(cfg.rutas.resultados, 'ajuste_prueba'); end
dirC = fullfile(dirA, 'cache');
dirF = fullfile(dirA, 'figuras');
for d = {dirA, dirC, dirF}, mkdir_si(d{1}); end
dirArff = cfg.ajuste.dir_arff;
if isempty(dirArff), dirArff = fullfile(cfg.rutas.resultados, 'arff', sprintf('d%d', cfg.ajuste.duracion)); end
I = readtable(fullfile(dirArff, 'indice.csv'), 'TextType', 'char');
sel = strcmp(I.condicion, cfg.ajuste.condicion) & strcmp(I.esquema, cfg.ajuste.esquema);
Ii = I(sel & strcmp(I.alcance, 'intra'), :);
rutas = fullfile(dirArff, Ii.archivo);
sesiones = Ii.id;
rutaGlobal = fullfile(dirArff, I.archivo{sel & strcmp(I.alcance, 'global')});
if cfg.modo_prueba
    rutas = rutas(1:min(2, end)); sesiones = sesiones(1:numel(rutas));
end
ns = numel(rutas);
nw = usar_paralelo(cfg);
registrar(cfg, 'Etapa 6: ajuste de clasificadores, %s %s, %d muestras, %s, %d sesiones | conjunto %s -> %s', ...
    cfg.ajuste.condicion, cfg.ajuste.esquema, cfg.ajuste.duracion, cfg.ajuste.config, ns, ...
    cfg.ajuste.conjunto, dirA);

%% 1) Diagnóstico de convergencia
archDiag = fullfile(dirA, 'diagnostico_convergencia.csv');
if ~isfile(archDiag)
    t0 = tic;
    partes = cell(ns, 1);
    parfor (s = 1:ns, nw)
        T = diagnostico_convergencia(rutas{s}, cfg);
        T.sesion = repmat(sesiones(s), height(T), 1);
        partes{s} = T;
    end
    writetable(vertcat(partes{:}), archDiag);
    registrar(cfg, '  diagnóstico listo (%.1f min).', toc(t0)/60);
end
Dg = readtable(archDiag, 'TextType', 'char');
graficar_diagnostico(Dg, dirF);

%% 2) Búsqueda en rejilla
modelos = rejillas_hiperparametros(cfg);
registrar(cfg, '  búsqueda: %d modelos x %d sesiones.', numel(modelos), ns);
t0 = tic;
parfor (s = 1:ns, nw)
    buscar_en_sesion(rutas{s}, modelos, cfg, fullfile(dirC, ['busqueda_' sesiones{s} '.mat']));
end
registrar(cfg, '  búsqueda lista (%.1f min).', toc(t0)/60);
B = cargar_cache(dirC, 'busqueda_', sesiones);
B = B(strcmp(B.estado, 'ok'), :);
writetable(B, fullfile(dirA, 'busqueda_resultados.csv'));
[g, ids] = findgroups(B.id);
Res = table(ids, splitapply(@(x) x(1), B.clf, g), splitapply(@(x) x(1), B.sel, g), ...
    splitapply(@numel, B.acc, g), splitapply(@mean, B.acc, g), splitapply(@std, B.acc, g), ...
    splitapply(@mean, B.kappa, g), splitapply(@mean, B.mae, g), ...
    'VariableNames', {'id','clf','sel','n_sesiones','acc_media','acc_de','kappa_media','mae_media'});
Res = Res(Res.n_sesiones == ns, :);               % solo modelos evaluados en todas las sesiones
Res = sortrows(Res, {'acc_media','kappa_media'}, {'descend','descend'});
writetable(Res, fullfile(dirA, 'busqueda_resumen.csv'));
graficar_busqueda(B, dirF);

%% 3) Modelos finales: línea base, mejores ajustados y ensamble
clfs = {'kNN','SVM','MLP','RF'};
base = [especificar_modelo('kNN', struct()), especificar_modelo('SVM', struct()), ...
        especificar_modelo('MLP', struct()), especificar_modelo('RF', struct())];
mejores = [];
for c = 1:numel(clfs)
    id = Res.id(find(strcmp(Res.clf, clfs{c}), 1));
    mejores = [mejores, modelos(strcmp({modelos.id}, id))]; %#ok<AGROW>
end
voto = especificar_modelo('Vote', struct('miembros', strjoin({mejores.clf}, '+')));
voto.miembros = mejores;
finales = [base, mejores, voto];
etiquetas = [strcat(clfs, ' (artículo)'), strcat(clfs, ' (ajustado)'), {'Ensamble (Vote)'}];
writetable(table(etiquetas(5:8).', {mejores.id}.', 'VariableNames', {'modelo','especificacion'}), ...
    fullfile(dirA, 'hiperparametros_elegidos.csv'));
registrar(cfg, '  mejores: %s', strjoin({mejores.id}, ' | '));

%% 4) Estimación final: CV repetida
semillas = cfg.ajuste.semillas_final;
t0 = tic;
parfor (s = 1:ns, nw)
    evaluar_finales(rutas{s}, finales, etiquetas, cfg, semillas, ...
        fullfile(dirC, ['final_' sesiones{s} '.mat']));
end
registrar(cfg, '  CV repetida lista (%.1f min).', toc(t0)/60);
F = cargar_cache(dirC, 'final_', sesiones);
F.p_binomial = arrayfun(@(a, n) p_binomial(round(a/100*n), n, 1/3), F.acc, F.n);
writetable(F, fullfile(dirA, 'final_por_sesion.csv'));

%% 5) Control de sesgo: CV anidada
ajust = F(contains(F.modelo, '(ajustado)'), :);
[g, mods] = findgroups(ajust.modelo);
[~, o] = sort(splitapply(@mean, ajust.acc, g), 'descend');
topA = mods(o(1:min(cfg.ajuste.n_anidada, end)));
An = table();
for t = 1:numel(topA)
    clf = strtok(topA{t});
    cand = rejillas_hiperparametros(cfg, clf);
    if numel(cand) > cfg.ajuste.anidada_max_config
        cand = cand(unique(round(linspace(1, numel(cand), cfg.ajuste.anidada_max_config))));
    end
    registrar(cfg, '  CV anidada de %s con %d candidatos internos.', clf, numel(cand));
    t0 = tic;
    filas = cell(ns, 1);
    parfor (s = 1:ns, nw)
        filas{s} = anidada_sesion(rutas{s}, cand, cfg, sesiones{s}, topA{t}, ...
            fullfile(dirC, sprintf('anidada_%s_%s.mat', clf, sesiones{s})));
    end
    An = [An; vertcat(filas{:})]; %#ok<AGROW>
    registrar(cfg, '    lista (%.1f min).', toc(t0)/60);
end
if ~isempty(An), writetable(An, fullfile(dirA, 'anidada_por_sesion.csv')); end

%% 6) Estadística y tabla para publicación
tabla = tabla_publicacion(F, An, etiquetas);
writetable(tabla, fullfile(dirA, 'tabla_ajuste_kmi3c.csv'));
graficar_final(F, etiquetas, dirF);
ganador = tabla.modelo{find(~contains(tabla.modelo, 'artículo'), 1)};
registrar(cfg, '  ganador (CV repetida): %s = %.2f %%', ganador, tabla.acc_media(strcmp(tabla.modelo, ganador)));

archPerm = fullfile(dirA, 'permutacion_global.csv');
if ~isfile(archPerm) && ~cfg.modo_prueba
    m = finales(strcmp(etiquetas, ganador));
    np = cfg.ajuste.n_permutaciones;
    obs = weka_validar(cargar_datos_ajuste(rutaGlobal, cfg), weka_clasificador_param(m, cfg), cfg).acc;
    perm = zeros(np, 1);
    t0 = tic;
    parfor (k = 1:np, nw)
        perm(k) = acc_permutada(rutaGlobal, m, cfg, k);
    end
    p = (1 + nnz(perm >= obs)) / (1 + np);
    writetable(table({ganador}, obs, mean(perm), std(perm), max(perm), np, p, 'VariableNames', ...
        {'modelo','acc_observada','acc_perm_media','acc_perm_de','acc_perm_max','n_perm','p'}), archPerm);
    registrar(cfg, '  permutación (global): acc = %.2f %%, nula = %.2f +- %.2f, p = %.3f (%.1f min).', ...
        obs, mean(perm), std(perm), p, toc(t0)/60);
end
registrar(cfg, 'Etapa 6 lista: %s', dirA);
end

%% ======================================================================
function mkdir_si(d)
if ~isfolder(d), mkdir(d); end
end

function T = cargar_cache(dirC, prefijo, sesiones)
partes = cell(numel(sesiones), 1);
for s = 1:numel(sesiones)
    a = fullfile(dirC, [prefijo sesiones{s} '.mat']);
    if ~isfile(a), continue, end
    S = load(a, 'R');
    R = S.R;
    R.sesion = repmat(sesiones(s), height(R), 1);
    partes{s} = R;
end
T = vertcat(partes{:});
end

function evaluar_finales(ruta, finales, etiquetas, cfg, semillas, ruta_cache)
if isfile(ruta_cache)                     % reutilizar solo si son los mismos modelos y semillas
    S = load(ruta_cache, 'R');
    if isequal(sort(S.R.especificacion), sort({finales.id}.')) && ...
            isequal(numel(S.R.acc_de_semillas), numel(finales)) && isfield(S, 'R') && ...
            all(strcmp(S.R.modelo, etiquetas(:))) && isequal(S.R.Properties.UserData, semillas)
        return
    end
end
D = cargar_datos_ajuste(ruta, cfg);
filas = cell(numel(finales), 1);
for i = 1:numel(finales)
    [Mm, Ms] = weka_cv_repetida(D, finales(i), cfg, semillas);
    filas{i} = table(etiquetas(i), {finales(i).id}, Mm.n, Mm.acc, Ms.acc, Mm.kappa, Mm.mae, ...
        Mm.wtpr, Mm.wfpr, Mm.wprec, Mm.wf, Mm.wroc, Mm.tpr, Mm.fpr, Mm.prec, {mat2str(Mm.confusion)}, ...
        'VariableNames', {'modelo','especificacion','n','acc','acc_de_semillas','kappa','mae', ...
        'wtpr','wfpr','wprec','wf','wroc','tpr','fpr','prec','confusion'});
end
R = vertcat(filas{:});
R.Properties.UserData = semillas;
tmp = [ruta_cache '.tmp.mat']; save(tmp, 'R'); movefile(tmp, ruta_cache, 'f');
end

function fila = anidada_sesion(ruta, cand, cfg, sesion, modelo, ruta_cache)
if isfile(ruta_cache)
    S = load(ruta_cache, 'fila'); fila = S.fila; return
end
D = cargar_datos_ajuste(ruta, cfg);
[M, eleg] = weka_cv_anidada(D, cand, cfg);
fila = table({sesion}, {modelo}, M.acc, M.kappa, M.mae, {strjoin(eleg, ' || ')}, ...
    numel(unique(eleg)), 'VariableNames', {'sesion','modelo','acc_anidada','kappa_anidada', ...
    'mae_anidada','elegidos_por_fold','n_distintos'});
tmp = [ruta_cache '.tmp.mat']; save(tmp, 'fila'); movefile(tmp, ruta_cache, 'f');
end

function a = acc_permutada(ruta, m, cfg, k)
D = cargar_datos_ajuste(ruta, cfg);
n = D.numInstances();
y = zeros(n, 1);
for i = 1:n, y(i) = D.instance(i-1).classValue(); end
s = RandStream('mt19937ar', 'Seed', 1000 + k);
y = y(randperm(s, n));
for i = 1:n, D.instance(i-1).setClassValue(y(i)); end
a = weka_validar(D, weka_clasificador_param(m, cfg), cfg).acc;
end

function p = p_binomial(k, n, p0)
% P(X >= k), X ~ Binomial(n, p0)  (betainc de MATLAB base)
if k <= 0, p = 1; else, p = betainc(p0, k, n - k + 1); end
end

function T = tabla_publicacion(F, An, etiquetas)
ses = unique(F.sesion);
filas = cell(numel(etiquetas), 1);
for i = 1:numel(etiquetas)
    Fi = F(strcmp(F.modelo, etiquetas{i}), :);
    a = Fi.acc;
    W = wilcoxon_holm([a, 100/3*ones(size(a))], {'m','azar'});
    p_base = NaN;
    if contains(etiquetas{i}, '(ajustado)')
        b = F(strcmp(F.modelo, strrep(etiquetas{i}, '(ajustado)', '(artículo)')), :);
        [~, ia, ib] = intersect(Fi.sesion, b.sesion);
        Wb = wilcoxon_holm([a(ia), b.acc(ib)], {'aj','base'});
        p_base = Wb.p;
    end
    acc_an = NaN; acc_an_de = NaN;
    if ~isempty(An) && any(strcmp(An.modelo, etiquetas{i}))
        x = An.acc_anidada(strcmp(An.modelo, etiquetas{i}));
        acc_an = mean(x); acc_an_de = std(x);
    end
    filas{i} = table(etiquetas(i), mean(a), std(a), min(a), max(a), mean(Fi.kappa), mean(Fi.mae), ...
        100*mean(Fi.wtpr), 100*mean(Fi.wfpr), 100*mean(Fi.wprec), 100*mean(Fi.wf), 100*mean(Fi.wroc), ...
        nnz(Fi.p_binomial < 0.05), numel(ses), W.p, p_base, acc_an, acc_an_de, ...
        'VariableNames', {'modelo','acc_media','acc_de','acc_min','acc_max','kappa','mae','TPR','FPR', ...
        'precision','F','ROC','sesiones_sobre_azar','n_sesiones','p_vs_azar','p_vs_base', ...
        'acc_anidada','acc_anidada_de'});
end
T = vertcat(filas{:});
aj = contains(T.modelo, '(ajustado)');          % Holm entre las 4 comparaciones vs base
p = T.p_vs_base(aj); [ps, o] = sort(p); adj = zeros(size(p)); acum = 0;
for j = 1:numel(ps), acum = max(acum, min(1, (numel(ps)-j+1)*ps(j))); adj(o(j)) = acum; end
T.p_vs_base_holm = nan(height(T), 1); T.p_vs_base_holm(aj) = adj;
T = sortrows(T, 'acc_media', 'descend');
end

function graficar_diagnostico(Dg, dirF)
f = figure('Visible', 'off', 'Position', [100 100 1200 380]);
M = Dg(strcmp(Dg.prueba, 'mlp_epocas'), :);
confs = unique(M.parametro, 'stable');
subplot(1,3,1); hold on
for i = 1:numel(confs)
    s = M(strcmp(M.parametro, confs{i}), :);
    [g, N] = findgroups(s.valor);
    plot(N, splitapply(@mean, s.acc_cv, g), '-o', 'DisplayName', ['CV ' confs{i}]);
    plot(N, splitapply(@mean, s.acc_entrenamiento, g), ':x', 'DisplayName', ['entren. ' confs{i}]);
end
set(gca, 'XScale', 'log'); grid on; xlabel('Épocas'); ylabel('Precisión (%)');
title('MLP: entrenamiento vs CV'); legend('FontSize', 6, 'Location', 'eastoutside');
subplot(1,3,2);
s = Dg(strcmp(Dg.prueba, 'rf_arboles'), :); [g, N] = findgroups(s.valor);
plot(N, splitapply(@mean, s.acc_cv, g), '-o'); grid on
xlabel('Árboles'); ylabel('Precisión CV (%)'); title('RF: n.º de árboles');
subplot(1,3,3);
s = Dg(strcmp(Dg.prueba, 'conjunto'), :);
[g, cj, cl] = findgroups(s.parametro, s.clasificador);
T = table(cj, cl, splitapply(@mean, s.acc_cv, g));
U = unstack(T, 'Var3', 'cl');
bar(categorical(U.cj), U{:, 2:end}); legend(U.Properties.VariableNames(2:end)); grid on
ylabel('Precisión CV (%)'); title('Familias de Hurst (parámetros del artículo)');
exportgraphics(f, fullfile(dirF, 'diagnostico_convergencia.png'), 'Resolution', 150);
close(f);
end

function graficar_busqueda(B, dirF)
f = figure('Visible', 'off', 'Position', [100 100 1200 380]);
k = B(strcmp(B.clf, 'kNN') & strcmp(B.sel, 'ninguna'), :);
k.K = cellfun(@(s) str2double(regexp(s, 'K=(\d+)', 'tokens', 'once')), k.id);
k.dist = cellfun(@(s) regexp(s, 'dist=(\w+)', 'tokens', 'once'), k.id);
k.peso = cellfun(@(s) regexp(s, 'peso=(\w+)', 'tokens', 'once'), k.id);
subplot(1,3,1); hold on
[g, d, w] = findgroups(k.dist, k.peso);
for i = 1:max(g)
    s = k(g == i, :); [gg, K] = findgroups(s.K);
    plot(K, splitapply(@mean, s.acc, gg), '-o', 'DisplayName', [d{i} '/' w{i}]);
end
grid on; xlabel('k'); ylabel('Precisión media (%)'); title('kNN (sin selección)'); legend('FontSize', 6);
v = B(strcmp(B.clf, 'SVM') & strcmp(B.sel, 'ninguna'), :);
v.C = cellfun(@(s) str2double(regexp(s, 'C=([\d.e-]+)', 'tokens', 'once')), v.id);
v.ker = cellfun(@(s) regexprep(regexprep(s, '^SVM\|', ''), ';?C=[\d.e-]+', ''), v.id, 'UniformOutput', false);
v.ker = regexprep(v.ker, '\|ninguna$', '');
subplot(1,3,2); hold on
[g, kk] = findgroups(v.ker);
for i = 1:max(g)
    s = v(g == i, :); [gg, C] = findgroups(s.C);
    plot(C, splitapply(@mean, s.acc, gg), '-o', 'DisplayName', kk{i});
end
set(gca, 'XScale', 'log'); grid on; xlabel('C'); ylabel('Precisión media (%)');
title('SVM (sin selección)'); legend('FontSize', 5, 'Location', 'eastoutside');
subplot(1,3,3);
[g, c, s] = findgroups(B.clf, B.sel);
T = table(c, s, accumarray_media(B, g));
U = unstack(T, 'Var3', 's');
bar(categorical(U.c), U{:, 2:end}); grid on; legend(U.Properties.VariableNames(2:end), 'FontSize', 6);
ylabel('Mejor precisión media (%)'); title('Selección de características');
exportgraphics(f, fullfile(dirF, 'busqueda_rejilla.png'), 'Resolution', 150);
close(f);
end

function v = accumarray_media(B, g)
% Para cada grupo (clf, sel): la mejor precisión media entre sesiones de sus modelos.
[gi, ~] = findgroups(B.id);
med = splitapply(@mean, B.acc, gi);
mejor = zeros(max(g), 1);
for k = 1:max(g)
    ids = unique(gi(g == k));
    mejor(k) = max(med(ids));
end
v = mejor;
end

function graficar_final(F, etiquetas, dirF)
ses = unique(F.sesion);
A = nan(numel(ses), numel(etiquetas));
for i = 1:numel(etiquetas)
    for s = 1:numel(ses)
        x = F.acc(strcmp(F.modelo, etiquetas{i}) & strcmp(F.sesion, ses{s}));
        if ~isempty(x), A(s, i) = x; end
    end
end
f = figure('Visible', 'off', 'Position', [100 100 1200 450]);
subplot(1,2,1);
grupo = categorical(repmat(etiquetas, numel(ses), 1), etiquetas);
boxchart(grupo(:), A(:)); hold on; yline(100/3, '--', 'azar'); grid on
ylabel('Precisión (CV 5-fold x 10)'); title('KMI 3 clases por sesión');
subplot(1,2,2);
[~, j] = max(mean(A, 1, 'omitnan'));
C = zeros(3);
Fi = F(strcmp(F.modelo, etiquetas{j}), :);
for i = 1:height(Fi)
    Ci = eval(Fi.confusion{i}); C = C + Ci ./ sum(Ci, 2);
end
heatmap({'N00','N10','N40'}, {'N00','N10','N40'}, round(100*C/height(Fi), 1), ...
    'XLabel', 'Predicha', 'YLabel', 'Real', 'Title', ['Confusión media: ' etiquetas{j}]);
exportgraphics(f, fullfile(dirF, 'final_por_sesion.png'), 'Resolution', 150);
close(f);
end
