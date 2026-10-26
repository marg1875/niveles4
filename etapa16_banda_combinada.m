function etapa16_banda_combinada(cfg)
%ETAPA16_BANDA_COMBINADA  Potencia de banda optimizada, combinada con escalamiento y entropía.
%   Con los 56 registros (14 sesiones KMI) la potencia de banda theta/mu/
%   beta y FuzzyEn quedaron arriba (sección 17 de docs/RESULTADOS.md), pero
%   la potencia de banda nunca se había optimizado. Aquí:
%   B1) Familias espectrales por canal, desde la PSD de Welch ya calculada
%       en las curvas de la etapa 13 (ventanas de 4 s, resolución ~1 Hz):
%         ABS  log potencia absoluta en delta, theta, mu, beta1, beta2, gamma
%         REL  potencia relativa (sobre 1-45 Hz) en las mismas bandas
%         RAZ  log cocientes theta/beta, mu/beta, theta/mu, (theta+mu)/beta, alfa/delta
%         FIN  potencia relativa en bandas finas de 2 Hz (4-30 Hz)
%         FOR  forma: frecuencia del pico alfa (7-13 Hz), borde espectral 90 %, centroide
%       y familias de escalamiento/entropía de las etapas 13-14 (FuzzyEn,
%       MFDFA h(+1), DFA3, wavelets, Higuchi).
%   B2) LOSO KMI 10 vs 40 (SVM ajustado ig8) de cada familia y de conjuntos
%       estándar (ABS/REL theta+mu+beta, todas las bandas, todos los cocientes).
%   B3) Selección hacia adelante entre familias (máx. cfg.banda.max_familias),
%       global y ANIDADA ESTRICTA: para cada sesión de prueba se repite toda
%       la selección con LOSO interna sobre las demás sesiones.
%   C)  Validación de finalistas: permutación por sesión de 1000 completas
%       (sin paro temprano), LOPO con permutación exacta por paciente (32),
%       IC 95 % por bootstrap de pacientes, voto por registro y decisión por
%       sesión. Finalistas: los mejores de B2/B3 y las configuraciones de la
%       sección 17 (potencia theta/mu/beta, FuzzyEn, MFDFA, DFA3, Hurst).
%   -> resultados/banda_combinada/*.csv

b = cfg.banda;
dirB = fullfile(cfg.rutas.resultados, b.carpeta);
dirC = fullfile(dirB, 'cache');
if ~isfolder(dirC), mkdir(dirC); end
nw = usar_paralelo(cfg);
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
m = modelo_desde_id(Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')});
registrar(cfg, '==== Etapa 16: potencia de banda optimizada y combinada (%d workers) ====', nw);
avance(dirB, 'pasos', 0, 4);

%% Matriz de familias (KMI, ventanas de 4 s)
[T, X, nombres, fam] = familias(cfg, b);
sesL = unique(T.sesion, 'stable');
C = numel(cfg.canales.nombres);
cols = @(conj) find(ismember(fam, strsplit(conj, '+')));
registrar(cfg, '  %d familias, %d atributos, %d ventanas KMI, %d sesiones', numel(unique(fam)), size(X, 2), height(T), numel(sesL));

%% B2) LOSO de cada familia y de conjuntos estándar
uf = unique(fam, 'stable');
estandar = {'ABSth+ABSmu+ABSb1+ABSb2', 'RELth+RELmu+RELb1+RELb2', strjoin(uf(startsWith(uf, 'ABS')), '+'), ...
    strjoin(uf(startsWith(uf, 'REL')), '+'), strjoin(uf(startsWith(uf, 'RAZ')), '+'), strjoin(uf(startsWith(uf, 'FIN')), '+'), ...
    strjoin(uf(startsWith(uf, 'FOR')), '+')};
conj = [uf, estandar];
L = loso_conjuntos(conj, T, X, nombres, cols, m, cfg, nw, dirC, dirB);
R = resumir(L);
writetable(R, fullfile(dirB, 'loso_familias.csv'));
registrar(cfg, '  B2: mejor familia/conjunto %s = %.1f %% (%d/%d sesiones > 50)', R.conjunto{1}, R.loso(1), R.sobre50(1), numel(sesL));
avance(dirB, 'pasos', 1, 4);

%% B3) Selección hacia adelante global y anidada estricta
pool = uf;
archA = fullfile(dirB, 'adelante.csv');
if ~isfile(archA)
    t0 = tic;
    res = cell(numel(sesL) + 1, 1);
    parfor (q = 1:numel(sesL) + 1, nw)
        entr = sesL; s = '';
        if q > 1, s = sesL{q - 1}; entr = sesL(~strcmp(sesL, s)); end
        res{q} = adelante(T, X, nombres, cols, pool, entr, s, m, cfg, b.max_familias);
    end
    A = vertcat(res{:});
    writetable(A, archA);
    registrar(cfg, '  B3: adelante %.1f min', toc(t0)/60);
end
A = readtable(archA, 'TextType', 'char');
g = A(strcmp(A.sesion_prueba, 'global'), :);
an = A(~strcmp(A.sesion_prueba, 'global'), :);
Wt = wilcoxon_holm([an.acc_prueba, 50 * ones(height(an), 1)], {'a', 'b'});
[u, ~, ic] = unique(an.conjunto); [nmax, im] = max(accumarray(ic, 1));
registrar(cfg, '  B3: global %s (LOSO %.1f %%) | ANIDADO ESTRICTO %.1f %% (%d/%d > 50, p_W = %.3f) | más elegido %s (%d/%d)', ...
    g.conjunto{1}, g.acc_interna(1), mean(an.acc_prueba), nnz(an.acc_prueba > 50), height(an), Wt.p, u{im}, nmax, height(an));
avance(dirB, 'pasos', 2, 4);

%% C) Validación completa de finalistas
fin = unique([R.conjunto(1:min(b.top_validar, height(R))).', g.conjunto(1), u(im), b.referencias], 'stable');
fin = fin(cellfun(@(c) all(ismember(strsplit(c, '+'), uf)), fin));
V = validar(fin, T, X, nombres, cols, m, cfg, nw, dirC, dirB, b);
writetable(V, fullfile(dirB, 'validacion_final.csv'));
avance(dirB, 'pasos', 4, 4);
registrar(cfg, '==== Etapa 16 lista: %s ====', dirB);
end

%% ======================================================================
function [T, X, nombres, fam] = familias(cfg, b)
% Familias espectrales (desde la PSD en caché) + escalamiento/entropía.
S = load(fullfile(cfg.rutas.resultados, b.dir_w, 'w_W1000.mat'));
W = S.W(strcmp(S.W.condicion, 'KMI'), :);
Cand = candidatos_rango(512, cfg);
baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for r = {'Re64a2', 'Re128a4'}
    c = Cand(strcmp({Cand.nombre}, r{1}));
    baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<AGROW>
end
[K, lay] = caract_registro(cfg, 1000, baseH);
M = double(W.feat);
ch = cfg.canales.nombres; C = numel(ch);
ib = strcmp({lay.bloques.nombre}, 'psd');
bl = lay.bloques(ib); f = bl.arg.f;
P = reshape(M(:, bl.o + (1:bl.K * C)), [], C, bl.K);           % ventanas x canal x frecuencia
pot = @(f1, f2) sum(P(:, :, f >= f1 & f < f2), 3);
tot = pot(1, 45);
bandas = {'de', 1, 4; 'th', 4, 8; 'mu', 8, 13; 'b1', 13, 20; 'b2', 20, 30; 'ga', 30, 45};
X = []; nombres = {}; fam = {};
    function agregar(nf, v)
        X = [X, v]; nombres = [nombres, strcat(nf, '_', ch)]; fam = [fam, repmat({nf}, 1, C)]; %#ok<SETNU>
    end
for k = 1:size(bandas, 1)
    pb = pot(bandas{k, 2}, bandas{k, 3});
    agregar(['ABS' bandas{k, 1}], log(pb));
    agregar(['REL' bandas{k, 1}], pb ./ tot);
end
th = pot(4, 8); mu = pot(8, 13); be = pot(13, 30); de = pot(1, 4);
agregar('RAZthb', log(th ./ be)); agregar('RAZmub', log(mu ./ be)); agregar('RAZthmu', log(th ./ mu));
agregar('RAZthmub', log((th + mu) ./ be)); agregar('RAZade', log(mu ./ de));
for f1 = 4:2:28
    agregar(sprintf('FIN%d', f1), pot(f1, f1 + 2) ./ tot);
end
ia = find(f >= 7 & f <= 13); [~, j] = max(P(:, :, ia), [], 3); agregar('FORiaf', f(ia(j)));
cs = cumsum(P(:, :, f >= 1 & f < 45), 3); fr = f(f >= 1 & f < 45);
[~, j] = max(cs >= 0.9 * cs(:, :, end), [], 3); agregar('FORsef', fr(j));
w = reshape(fr, 1, 1, []); agregar('FORcen', sum(P(:, :, f >= 1 & f < 45) .* w, 3) ./ tot);
% escalamiento y entropía (nombres de la etapa 13, sin el prefijo W1000)
[~, ik] = ismember(strcat('W1000', b.fractales), {K.nombre});
Xf = caract_desde_curvas(M, K(ik), lay, cfg);
for k = 1:numel(b.fractales)
    agregar(b.fractales{k}, Xf(:, (k-1)*C + (1:C)));
end
X(~isfinite(X)) = NaN;
T = removevars(W, intersect({'feat', 'valida'}, W.Properties.VariableNames));
end

function Lt = loso_conjuntos(conj, T, X, nombres, cols, m, cfg, nw, dirC, dirB)
arch = fullfile(dirC, strcat('L_', strrep(conj, '+', '__'), '.mat'));
pend = find(~cellfun(@isfile, arch));
if ~isempty(pend)
    t0 = tic;
    Xg = cellfun(@(c) X(:, cols(c)), conj(pend), 'UniformOutput', false);
    Ng = cellfun(@(c) nombres(cols(c)), conj(pend), 'UniformOutput', false);
    ap = arch(pend); cp = conj(pend);
    parfor (q = 1:numel(pend), nw)
        Tq = T; Tq.feat = Xg{q};
        acc = loso_por_sesion(Tq, Ng{q}, m, cfg);
        R = table(repmat(cp(q), numel(acc), 1), unique(T.sesion, 'stable'), acc, 'VariableNames', {'conjunto','sesion','loso'});
        guardar(ap{q}, R);
    end
    registrar(cfg, '  LOSO de %d conjuntos (%.1f min)', numel(pend), toc(t0)/60);
end
Lt = cellfun(@(a) getfield(load(a), 'R'), arch, 'UniformOutput', false); %#ok<GFLD>
Lt = vertcat(Lt{:});
avance(dirB, 'loso', numel(conj), numel(conj));
end

function R = resumir(L)
[g, cj] = findgroups(L.conjunto);
R = table(cj, splitapply(@(x) mean(x, 'omitnan'), L.loso, g), splitapply(@(x) nnz(x > 50), L.loso, g), ...
    'VariableNames', {'conjunto','loso','sobre50'});
R = sortrows(R, 'loso', 'descend');
end

function F = adelante(T, X, nombres, cols, pool, entr, sprueba, m, cfg, maxf)
% Selección hacia adelante con LOSO interna en 'entr'; si sprueba no está
% vacía, se evalúa el conjunto elegido en esa sesión (entrenando sin ella).
weka_iniciar(cfg);
k = ismember(T.sesion, entr);
interna = @(c) mean(loso_por_sesion(sub(T(k, :), X(k, cols(c))), nombres(cols(c)), m, cfg), 'omitnan');
v = cellfun(interna, pool);
[best, j] = max(v); actual = pool(j); restantes = pool; restantes(j) = [];
while numel(actual) < maxf && ~isempty(restantes)
    cand = cellfun(@(r) strjoin([actual, {r}], '+'), restantes, 'UniformOutput', false);
    v = cellfun(interna, cand);
    [vm, j] = max(v);
    if vm <= best, break, end
    best = vm; actual = [actual, restantes(j)]; restantes(j) = []; %#ok<AGROW>
end
conj = strjoin(actual, '+');
acc = NaN; et = 'global';
if ~isempty(sprueba)
    et = sprueba;
    a = loso_por_sesion(sub(T, X(:, cols(conj))), nombres(cols(conj)), m, cfg, {sprueba});
    acc = a(strcmp(unique(T.sesion, 'stable'), sprueba));
end
F = table({et}, {conj}, best, acc, 'VariableNames', {'sesion_prueba','conjunto','acc_interna','acc_prueba'});
end

function T = sub(T, X)
T.feat = X;
end

function V = validar(fin, T, X, nombres, cols, m, cfg, nw, dirC, dirB, b)
% Permutación completa (1000) por sesión, LOPO con permutación exacta por
% paciente, IC por bootstrap de pacientes, voto por registro.
pac = regexprep(T.sesion, '_.*$', '');
filas = cell(numel(fin), 1);
for i = 1:numel(fin)
    arch = fullfile(dirC, sprintf('val_%02d.mat', i));
    if isfile(arch), S = load(arch); filas{i} = S.fila; continue, end
    t0 = tic;
    Tk = sub(T, X(:, cols(fin{i}))); nom = nombres(cols(fin{i}));
    pS = predicciones_grupo(Tk, nom, m, cfg, Tk.sesion);
    pP = predicciones_grupo(Tk, nom, m, cfg, pac);
    [aS, ses] = acc_sesion(Tk, pS); aP = acc_sesion(Tk, pP);
    pacS = regexprep(ses, '_.*$', '');
    % permutación por sesión: 1000 completas
    sesC = unique(Tk.sesion); nula = zeros(b.n_perm, 1);
    for b0 = 1:6*max(nw,1):b.n_perm
        idx = b0:min(b.n_perm, b0 + 6*max(nw,1) - 1); parte = zeros(numel(idx), 1);
        parfor (j = 1:numel(idx), nw)
            rs = RandStream('mt19937ar', 'Seed', 500 + idx(j));
            vol = sesC(rand(rs, numel(sesC), 1) < 0.5);
            parte(j) = mean(acc_sesion(voltear(Tk, vol), predicciones_grupo(Tk, nom, m, cfg, Tk.sesion, vol)));
        end
        nula(idx) = parte;
        avance(dirB, sprintf('perm_%02d', i), idx(end), b.n_perm);
    end
    pSes = (1 + nnz(nula >= mean(aS))) / (1 + b.n_perm);
    % permutación exacta por paciente (LOPO)
    up = unique(pac, 'stable'); nulaP = zeros(2^numel(up), 1);
    parfor (q = 1:2^numel(up), nw)
        vol = unique(Tk.sesion(ismember(pac, up(bitget(q - 1, 1:numel(up)) == 1))));
        nulaP(q) = mean(acc_sesion(voltear(Tk, vol), predicciones_grupo(Tk, nom, m, cfg, pac, vol)));
    end
    pPac = mean(nulaP >= mean(aP));
    [vA, vS] = votar(Tk, pS);
    ciS = boot(aS, pacS); ciP = boot(aP, pacS);
    fila = table(fin(i), numel(nom), mean(aS), ciS(1), ciS(2), nnz(aS > 50), numel(aS), pSes, mean(nula), ...
        mean(aP), ciP(1), ciP(2), pPac, vA(1), vA(2), vS(1), vS(2), ...
        'VariableNames', {'conjunto','n_atributos','loso','loso_ic_inf','loso_ic_sup','sobre50','n_ses', ...
        'p_perm_sesion_1000','nula_media','lopo','lopo_ic_inf','lopo_ic_sup','p_perm_paciente', ...
        'voto_archivo','voto_archivo_p','decision_sesion','decision_sesion_p'});
    save(arch, 'fila'); filas{i} = fila;
    registrar(cfg, '  final %-40s LOSO %.1f %% [%.0f-%.0f] p1000 = %.3f | LOPO %.1f %% p_pac = %.3f | voto %.0f %% | sesión %.0f %% (p = %.3f) (%.1f min)', ...
        fin{i}, mean(aS), ciS, pSes, mean(aP), pPac, vA(1), vS(1), vS(2), toc(t0)/60);
end
V = sortrows(vertcat(filas{:}), 'loso', 'descend');
end

function T = voltear(T, vol)
v = ismember(T.sesion, vol); T.nivel(v) = 50 - T.nivel(v);
end

function [a, ses] = acc_sesion(T, pred)
ses = unique(T.sesion, 'stable'); a = nan(numel(ses), 1);
for s = 1:numel(ses)
    k = strcmp(T.sesion, ses{s}) & isfinite(pred);
    if numel(unique(T.nivel(k))) == 2, a(s) = 100 * mean(pred(k) == T.nivel(k)); end
end
ok = isfinite(a); a = a(ok); ses = ses(ok);
end

function ci = boot(a, pac)
rs = RandStream('mt19937ar', 'Seed', 7); up = unique(pac); mm = zeros(10000, 1);
for i = 1:10000
    sel = up(randi(rs, numel(up), numel(up), 1));
    mm(i) = mean(cell2mat(cellfun(@(p) a(strcmp(pac, p)), sel, 'UniformOutput', false)));
end
ci = quantile(mm, [0.025 0.975]);
end

function [A, S] = votar(T, pred)
[g, ses, ~] = findgroups(T.sesion, T.archivo);
f40 = splitapply(@(p) mean(p(isfinite(p)) == 40), pred, g);
niv = splitapply(@(x) x(1), T.nivel, g);
ac = (f40 > 0.5) == (niv == 40);
A = [100 * mean(ac), 1 - binocdf(nnz(ac) - 1, numel(ac), 0.5)];
us = unique(ses, 'stable'); d = nan(numel(us), 1);
for s = 1:numel(us)
    i10 = strcmp(ses, us{s}) & niv == 10; i40 = strcmp(ses, us{s}) & niv == 40;
    if any(i10) && any(i40), d(s) = f40(i40) > f40(i10); end
end
d = d(isfinite(d));
S = [100 * mean(d), 1 - binocdf(nnz(d) - 1, numel(d), 0.5)];
end

function guardar(ruta, R) %#ok<INUSD>
tmp = [ruta '.tmp.mat']; save(tmp, 'R'); movefile(tmp, ruta, 'f');
end

function avance(dirB, nombre, hecho, total)
fid = fopen(fullfile(dirB, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end
