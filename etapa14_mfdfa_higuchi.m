function etapa14_mfdfa_higuchi(cfg)
%ETAPA14_MFDFA_HIGUCHI  MFDFA h(+1) combinado con Higuchi: configuración ganadora.
%   Ventanas de cfg.mfh.L muestras (las curvas de la etapa 13). Pasos:
%   A) LOSO KMI 10 vs 40 (SVM ig8), por sesión, de TODAS las variantes de
%      las familias de cfg.mfh.familias (MFDFA h(+1), Higuchi y controles).
%   B) Pares MFDFA h(+1) x Higuchi (las cfg.mfh.top mejores de cada una por
%      LOSO), pares de MFDFA con los controles y tríos MF + HIG + control.
%   C) LOSO anidada ESTRICTA de todo el procedimiento: para cada sesión de
%      prueba s, con las otras 12 sesiones se calcula la LOSO interna de
%      todas las variantes de MF y HIG, se toman las cfg.mfh.top mejores de
%      cada una, se evalúan sus pares, se elige el mejor (individual o par)
%      y se mide en s. También las reglas "solo MF" y "solo HIG".
%   D) Validación final (evaluar_generalizacion, KMI y FIS) y permutación
%      secuencial de los mejores pares, las mejores individuales y la
%      configuración más elegida por el anidado; con los clasificadores
%      SVM ig8, ig16 y sin selección para los dos mejores pares.
%   -> resultados/mfdfa_higuchi/*.csv

p = cfg.mfh;
dirM = fullfile(cfg.rutas.resultados, 'mfdfa_higuchi');
dirC = fullfile(dirM, 'cache');
for d = {dirM, dirC}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
nw = usar_paralelo(cfg);
Cand = candidatos_rango(512, cfg);
baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for b = {'Re64a2', 'Re128a4'}
    c = Cand(strcmp({Cand.nombre}, b{1}));
    baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<AGROW>
end
L = p.L; pre = sprintf('W%d', L);
[K, lay] = caract_registro(cfg, L, baseH);
archW = fullfile(cfg.rutas.resultados, 'complejidad', 'cache', ['w_' pre '.mat']);
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
spec = Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')};
m8 = modelo_desde_id(spec);
registrar(cfg, '==== Etapa 14: MFDFA h(+1) + Higuchi en %s (%d workers) ====', pre, nw);
avance(dirM, 'pasos', 0, 4);

%% Matriz de todas las variantes de las familias de interés (KMI)
fam = {K.familia};
sel = ismember(fam, p.familias);
Ks = K(sel);
S = load(archW); WK = S.W(strcmp(S.W.condicion, 'KMI'), :);
[X, nombres] = caract_desde_curvas(double(WK.feat), Ks, lay, cfg);
meta = removevars(WK, 'feat');
sesL = unique(meta.sesion, 'stable');
C = lay.C;
cols = @(conj) cell2mat(cellfun(@(c) (find(strcmp({Ks.nombre}, c)) - 1) * C + (1:C), strsplit(conj, '+'), 'UniformOutput', false));
nomDe = @(conj) nombres(cols(conj));
registrar(cfg, '  %d variantes (%s), %d ventanas KMI', numel(Ks), strjoin(p.familias, ', '), height(meta));

%% A) LOSO por sesión de todas las variantes
indiv = {Ks.nombre};
Lind = loso_conjuntos(indiv, X, meta, cols, nomDe, m8, cfg, nw, dirC, dirM, 'A');
writetable(resumir(Lind, Ks), fullfile(dirM, 'individuales.csv'));
avance(dirM, 'pasos', 1, 4);

%% B) Pares y tríos
famDe = containers.Map({Ks.nombre}, {Ks.familia});
topDe = @(f, n) top_por(Lind, indiv(strcmp(cellfun(@(c) famDe(c), indiv, 'UniformOutput', false), f)), n);
mf = topDe('MFhp1', p.top); hg = topDe('HIG', p.top);
combos = {};
for a = 1:numel(mf), for b = 1:numel(hg), combos{end+1} = [mf{a} '+' hg{b}]; end, end %#ok<AGROW>
for f = setdiff(p.familias, {'MFhp1', 'HIG'})
    ct = topDe(f{1}, 3);
    for a = 1:3, for b = 1:numel(ct), combos{end+1} = [mf{a} '+' ct{b}]; end, end %#ok<AGROW>
    combos{end+1} = [mf{1} '+' hg{1} '+' ct{1}]; %#ok<AGROW>
end
Lcom = loso_conjuntos(combos, X, meta, cols, nomDe, m8, cfg, nw, dirC, dirM, 'B');
Rc = resumir(Lcom, Ks);
writetable(Rc, fullfile(dirM, 'combinaciones.csv'));
registrar(cfg, '  mejor individual MF %s = %.1f %% | HIG %s = %.1f %% | mejor combinación %s = %.1f %%', ...
    mf{1}, media(Lind, mf{1}), hg{1}, media(Lind, hg{1}), Rc.conjunto{1}, Rc.loso(1));
avance(dirM, 'pasos', 2, 4);

%% C) LOSO anidada estricta del procedimiento completo
archN = fullfile(dirM, 'anidado.csv');
if ~isfile(archN)
    t0 = tic;
    res = cell(numel(sesL), 1);
    mfAll = indiv(strcmp(cellfun(@(c) famDe(c), indiv, 'UniformOutput', false), 'MFhp1'));
    hgAll = indiv(strcmp(cellfun(@(c) famDe(c), indiv, 'UniformOutput', false), 'HIG'));
    top = p.top;
    parfor (s = 1:numel(sesL), nw)
        res{s} = anidado_sesion(s, sesL, X, meta, cols, nomDe, mfAll, hgAll, top, m8, cfg);
    end
    N = vertcat(res{:});
    writetable(N, archN);
    for regla = {'MF+HIG', 'solo MF', 'solo HIG'}
        a = N.acc(strcmp(N.regla, regla{1}));
        Wt = wilcoxon_holm([a, 50 * ones(numel(a), 1)], {'a', 'b'});
        e = N.elegido(strcmp(N.regla, regla{1}));
        [u, ~, ic] = unique(e); [~, im] = max(accumarray(ic, 1));
        registrar(cfg, '  anidado estricto %-8s: LOSO %.1f %% (%d/%d > 50, p_W = %.3f) | más elegido %s (%d/%d)', ...
            regla{1}, mean(a), nnz(a > 50), numel(a), Wt.p, u{im}, max(accumarray(ic, 1)), numel(a));
    end
    registrar(cfg, '  anidado: %.1f min', toc(t0)/60);
end
N = readtable(archN, 'TextType', 'char');
avance(dirM, 'pasos', 3, 4);

%% D) Validación final y permutación secuencial
e = N.elegido(strcmp(N.regla, 'MF+HIG'));
[u, ~, ic] = unique(e); [~, im] = max(accumarray(ic, 1));
fin = unique([Rc.conjunto(1:min(5, height(Rc))).', mf(1), hg(1), u(im)], 'stable');
mods = {m8};
etq = {'SVM ig8'};
for s2 = {'ig16', 'ninguna'}
    mods{end+1} = modelo_desde_id(spec, s2{1}); etq{end+1} = ['SVM ' s2{1}]; %#ok<AGROW>
end
trab = {};
for i = 1:numel(fin), trab(end+1, :) = {fin{i}, 1}; end %#ok<AGROW>
for i = 1:min(2, height(Rc)), for k = 2:3, trab(end+1, :) = {Rc.conjunto{i}, k}; end, end %#ok<AGROW>
Wall = removevars(S.W, 'feat');
[Xall, ~] = caract_desde_curvas(double(S.W.feat), Ks, lay, cfg);
filas = cell(size(trab, 1), 1);
for i = 1:size(trab, 1)
    arch = fullfile(dirC, sprintf('final_%02d.mat', i));
    if isfile(arch), Sf = load(arch); filas{i} = Sf.fila; continue, end
    t0 = tic;
    cj = trab{i, 1}; mm = mods{trab{i, 2}};
    T = Wall; T.feat = Xall(:, cols(cj));
    G = evaluar_generalizacion(T, nomDe(cj), mm, cfg);
    TK = T(strcmp(T.condicion, 'KMI'), :);
    [obs, nula, pp] = permutacion_secuencial(TK, nomDe(cj), mm, cfg, nw, p.n_perm, 20, ...
        fullfile(dirM, sprintf('perm_%02d.avance.txt', i)));
    fila = [table({cj}, etq(trab{i, 2}), numel(cols(cj)), obs, mean(nula), numel(nula), pp, ...
        'VariableNames', {'conjunto','modelo','n_atributos','loso_perm','nula_media','n_perm','p_perm'}), G];
    save(arch, 'fila'); filas{i} = fila;
    registrar(cfg, '  final %-45s %-12s LOSO %.1f %% (p = %.4f, %d perm.) | %d/%d > 50 | FIS %.1f %% | Px010 %.0f %% (%.1f min)', ...
        cj, etq{trab{i, 2}}, obs, pp, numel(nula), G.ses_sobre50_kmi, G.n_ses_kmi, G.loso_fis, G.px010_invertido, toc(t0)/60);
end
V = sortrows(vertcat(filas{:}), 'loso_perm', 'descend');
writetable(V, fullfile(dirM, 'validacion_final.csv'));
avance(dirM, 'pasos', 4, 4);
registrar(cfg, '==== Etapa 14 lista: %s ====', dirM);
end

%% ======================================================================
function Lt = loso_conjuntos(conj, X, meta, cols, nomDe, m, cfg, nw, dirC, dirM, etiqueta)
% LOSO por sesión de cada conjunto, con caché por conjunto.
arch = fullfile(dirC, strcat('L_', strrep(conj, '+', '__'), '.mat'));
pend = find(~cellfun(@isfile, arch));
t0 = tic;
lote = 6 * max(nw, 1);
for b0 = 1:lote:numel(pend)
    g = pend(b0:min(numel(pend), b0 + lote - 1));
    Xg = cellfun(@(c) X(:, cols(c)), conj(g), 'UniformOutput', false);
    Ng = cellfun(@(c) nomDe(c), conj(g), 'UniformOutput', false);
    parfor (q = 1:numel(g), nw)
        T = meta; T.feat = Xg{q};
        acc = loso_por_sesion(T, Ng{q}, m, cfg);
        R = table(repmat(conj(g(q)), numel(acc), 1), unique(meta.sesion, 'stable'), acc, ...
            'VariableNames', {'conjunto','sesion','loso'});
        guardar(arch{g(q)}, R);
    end
    avance(dirM, ['loso_' etiqueta], b0 + numel(g) - 1, numel(pend));
end
if ~isempty(pend), registrar(cfg, '  LOSO %s: %d conjuntos (%.1f min)', etiqueta, numel(pend), toc(t0)/60); end
Lt = cellfun(@(a) getfield(load(a), 'R'), arch, 'UniformOutput', false); %#ok<GFLD>
Lt = vertcat(Lt{:});
end

function R = anidado_sesion(s, sesL, X, meta, cols, nomDe, mfAll, hgAll, top, m, cfg)
% Procedimiento completo sin la sesión s, evaluado en s.
weka_iniciar(cfg);
tr = ~strcmp(meta.sesion, sesL{s});
interna = @(c) mean(loso_por_sesion(sub(meta(tr, :), X(tr, cols(c))), nomDe(c), m, cfg), 'omitnan');
vMF = cellfun(interna, mfAll); vHG = cellfun(interna, hgAll);
[~, o1] = sort(vMF, 'descend'); [~, o2] = sort(vHG, 'descend');
mf = mfAll(o1(1:min(top, end))); hg = hgAll(o2(1:min(top, end)));
pares = {};
for a = 1:numel(mf), for b = 1:numel(hg), pares{end+1} = [mf{a} '+' hg{b}]; end, end %#ok<AGROW>
vP = cellfun(interna, pares);
cand = [mfAll, hgAll, pares]; v = [vMF, vHG, vP];
[~, j] = max(v); eleg = {cand{j}, mf{1}, hg{1}};
reglas = {'MF+HIG', 'solo MF', 'solo HIG'};
acc = zeros(3, 1);
for r = 1:3
    a = loso_por_sesion(sub(meta, X(:, cols(eleg{r}))), nomDe(eleg{r}), m, cfg, sesL(s));
    acc(r) = a(s);
end
R = table(repmat(sesL(s), 3, 1), reglas.', eleg.', acc, 'VariableNames', {'sesion','regla','elegido','acc'});
end

function T = sub(meta, X)
T = meta; T.feat = X;
end

function t = top_por(Lt, conj, n)
v = cellfun(@(c) media(Lt, c), conj);
[~, o] = sort(v, 'descend');
t = conj(o(1:min(n, numel(o))));
end

function v = media(Lt, c)
v = mean(Lt.loso(strcmp(Lt.conjunto, c)), 'omitnan');
end

function R = resumir(Lt, Ks)
[g, cj] = findgroups(Lt.conjunto);
R = table(cj, splitapply(@(x) mean(x, 'omitnan'), Lt.loso, g), splitapply(@(x) nnz(x > 50), Lt.loso, g), ...
    'VariableNames', {'conjunto','loso','sobre50'});
fam = cell(height(R), 1);
for i = 1:height(R)
    partes = strsplit(R.conjunto{i}, '+');
    [~, ik] = ismember(partes, {Ks.nombre});
    fam{i} = strjoin({Ks(ik).familia}, '+');
end
R.familias = fam;
R = sortrows(R, 'loso', 'descend');
end

function guardar(ruta, R) %#ok<INUSD>
tmp = [ruta '.tmp.mat']; save(tmp, 'R'); movefile(tmp, ruta, 'f');
end

function avance(dirM, nombre, hecho, total)
fid = fopen(fullfile(dirM, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end
