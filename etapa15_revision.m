function etapa15_revision(cfg)
%ETAPA15_REVISION  Controles pedidos por un revisor (puntos 3, 4, 5 y 7).
%   Configuraciones: las mejores de las etapas 13-14 (MFDFA h(+1), MFDFA +
%   Higuchi, DFA3), la referencia de Hurst, FuzzyEn (detector de huella) y
%   dos baselines estándar de BCI calculados aquí: CSP (8-30 Hz) +
%   regresión logística y potencia de banda (theta, mu, beta) por canal.
%   Ventanas contiguas de 4 s (las de la etapa 13; Hurst en 2 s).
%
%   3) Orden temporal frente a nivel: cada día (paciente x mes) tiene 4
%      registros (FIS10, KMI10, FIS40, KMI40) con su hora real. Se entrena
%      el clasificador 10/40 con la KMI de los OTROS días y se predicen los
%      4 registros del día. Si aprendiera "más tarde en el día", la
%      fracción de ventanas predichas 40 seguiría la hora (incluida la FIS
%      y el KMI10 frente al FIS10); si aprende el nivel, seguiría el nivel.
%      Modelo mixto: frac40 ~ nivel40 + minuto + FIS + (1 | día).
%   4) Validación dejando un PACIENTE fuera (LOPO), permutación exacta por
%      paciente (2^5 = 32 asignaciones; p mínimo 1/32) e intervalos de
%      confianza por bootstrap de pacientes (LOSO y LOPO).
%   5) Decisión por registro: voto mayoritario de las ventanas de cada
%      archivo (26 archivos KMI) y decisión por sesión ("¿cuál de los dos
%      registros es el de 40 %?", 13 sesiones), con prueba binomial exacta.
%   7) Baselines con los mismos controles (LOSO con permutación secuencial,
%      LOPO, voto y orden frente a nivel).
%   -> resultados/revision/*.csv

rv = cfg.revision;
dirR = fullfile(cfg.rutas.resultados, rv.carpeta);
dirC = fullfile(dirR, 'cache');
for d = {dirR, dirC}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
nw = usar_paralelo(cfg);
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
mSVM = modelo_desde_id(Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')});
mLOG = especificar_modelo('LOG', struct('ridge', 0.5));
registrar(cfg, '==== Etapa 15: controles de revisión (%d workers) ====', nw);
avance(dirR, 'pasos', 0, 5);

%% Configuraciones y sus matrices de características (todas las ventanas, KMI y FIS)
conf = {
    'MFDFA h(+1) o3 n5-24',      'W1000MFhp1o3n5a24',               'svm'
    'MFDFA h(+1) o3 n5-32',      'W1000MFhp1o3n5a32',               'svm'
    'MFDFA h(+1) + Higuchi',     'W1000MFhp1o3n6a24+W1000HIG12a32', 'svm'
    'DFA3 traslape n5-32',       'W1000DFAto5a32',                  'svm'
    'Hurst Re64a2+Re128a4 (2 s)','W512Re64a2+W512Re128a4',          'svm'
    'FuzzyEn r0.10 t2',          'W1000FEr10t2',                    'svm'
    'Potencia de banda (th/mu/b)','BANDAS',                          'svm'
    'Potencia de banda + LogReg','BANDAS',                          'log'
    'CSP 8-30 Hz + LogReg',      'CSP',                             'log'};
D = datos(conf, cfg, dirC, nw);
modelo = @(k) ifelse(strcmp(conf{k, 3}, 'svm'), mSVM, mLOG);
avance(dirR, 'pasos', 1, 5);

%% 4 y 5) LOSO y LOPO (predicciones por ventana), voto por registro, bootstrap y permutación por paciente
filas = {}; votos = {};
for k = 1:size(conf, 1)
    archK = fullfile(dirC, sprintf('val_%02d.mat', k));
    if isfile(archK)
        S = load(archK); filas(end+1, :) = S.fila; votos = [votos; S.voto]; continue %#ok<AGROW>
    end
    t0 = tic;
    T = D.(campo(conf{k, 2})); TK = T(strcmp(T.condicion, 'KMI'), :);
    [nom, ff] = config_feat(conf{k, 2}, D, TK, cfg);
    pac = regexprep(TK.sesion, '_.*$', '');
    pLOSO = predicciones_grupo(TK, nom, modelo(k), cfg, TK.sesion, {}, [], ff);
    pLOPO = predicciones_grupo(TK, nom, modelo(k), cfg, pac, {}, [], ff);
    [aS, sesS] = acc_sesion(TK, pLOSO); aP = acc_sesion(TK, pLOPO);
    pacS = regexprep(sesS, '_.*$', '');
    ciS = bootstrap_pac(aS, pacS); ciP = bootstrap_pac(aP, pacS);
    % permutación exacta por paciente (LOPO)
    up = unique(pac, 'stable'); nula = zeros(2^numel(up), 1);
    parfor (q = 1:2^numel(up), nw)
        vol = up(bitget(q - 1, 1:numel(up)) == 1);
        vs = unique(TK.sesion(ismember(pac, vol)));
        pq = predicciones_grupo(TK, nom, modelo(k), cfg, pac, vs, [], ff);
        nula(q) = mean(acc_sesion_vol(TK, pq, vs));
    end
    pPac = mean(nula >= mean(aP));
    % LOSO con permutación secuencial por sesión (solo baselines; el resto ya está en las etapas 13-14)
    pSes = NaN; nPerm = 0;
    if rv.perm_todas || ismember(conf{k, 2}, {'BANDAS', 'CSP'})
        [pSes, nPerm] = perm_secuencial(TK, nom, modelo(k), cfg, ff, mean(aS), nw, ...
            fullfile(dirR, sprintf('perm_%02d.avance.txt', k)));
    end
    % voto por registro y decisión por sesión
    [vA, vS, voto] = votar(TK, pLOSO, conf{k, 1}, 'LOSO');
    [vAP, vSP, votoP] = votar(TK, pLOPO, conf{k, 1}, 'LOPO');
    fila = {conf{k, 1}, conf{k, 2}, conf{k, 3}, mean(aS), ciS(1), ciS(2), pSes, nPerm, nnz(aS > 50), numel(aS), ...
        mean(aP), ciP(1), ciP(2), pPac, nnz(aP > 50), vA.acc, vA.p, vS.acc, vS.p, vAP.acc, vAP.p, vSP.acc, vSP.p};
    voto = [voto; votoP]; %#ok<AGROW>
    save(archK, 'fila', 'voto');
    filas(end+1, :) = fila; votos = [votos; voto]; %#ok<AGROW>
    registrar(cfg, ['  %-28s LOSO %.1f %% [%.1f-%.1f] | LOPO %.1f %% [%.1f-%.1f] p_pac = %.3f | ' ...
        'voto archivo %.0f %% (p = %.3f) | sesión %.0f %% (p = %.3f) (%.1f min)'], conf{k, 1}, mean(aS), ciS, ...
        mean(aP), ciP, pPac, vA.acc, vA.p, vS.acc, vS.p, toc(t0)/60);
end
V = cell2table(filas, 'VariableNames', {'configuracion','conjunto','clasificador','loso','loso_ic_inf','loso_ic_sup', ...
    'loso_p_perm_sesion','loso_n_perm','loso_sobre50','n_sesiones','lopo','lopo_ic_inf','lopo_ic_sup','lopo_p_perm_paciente', ...
    'lopo_sobre50','voto_archivo_loso','voto_archivo_loso_p','decision_sesion_loso','decision_sesion_loso_p', ...
    'voto_archivo_lopo','voto_archivo_lopo_p','decision_sesion_lopo','decision_sesion_lopo_p'});
writetable(V, fullfile(dirR, 'loso_lopo_voto.csv'));
writetable(votos, fullfile(dirR, 'votos_por_registro.csv'));
avance(dirR, 'pasos', 3, 5);

%% 3) Orden temporal frente a nivel (entrenar con KMI de otros días, predecir los 4 registros del día)
O = table(); M = {};
for k = 1:size(conf, 1)
    archO = fullfile(dirC, sprintf('orden_%02d.mat', k));
    if isfile(archO), S = load(archO); O = [O; S.Ok]; M(end+1, :) = S.fila; continue, end %#ok<AGROW>
    T = D.(campo(conf{k, 2}));
    dia = regexprep(T.sesion, '_(KMI|FIS)$', '');
    completos = dias_completos(T, dia);
    completos = setdiff(completos, rv.dias_excluidos_orden);   % legados: su hora de guardado no es la de grabación
    T = T(ismember(dia, completos), :); dia = dia(ismember(dia, completos));
    [nom, ff] = config_feat(conf{k, 2}, D, T, cfg);
    pr = predicciones_grupo(T, nom, modelo(k), cfg, dia, {}, strcmp(T.condicion, 'KMI'), ff);
    Ok = fraccion_por_registro(T, pr, dia, D.horas, conf{k, 1});
    [fila, txt] = modelo_mixto(Ok, conf{k, 1});
    save(archO, 'Ok', 'fila');
    O = [O; Ok]; M(end+1, :) = fila; %#ok<AGROW>
    registrar(cfg, '  orden vs nivel %-28s %s', conf{k, 1}, txt);
end
writetable(O, fullfile(dirR, 'orden_vs_nivel_por_registro.csv'));
writetable(cell2table(M, 'VariableNames', {'configuracion','b_nivel40','p_nivel40','b_minuto10','p_minuto10','b_fis','p_fis', ...
    'c_KMI10_menos_FIS10','c_KMI40_menos_FIS40','c_FIS40_menos_FIS10','c_KMI40_menos_KMI10','n_dias'}), ...
    fullfile(dirR, 'orden_vs_nivel_modelo.csv'));
avance(dirR, 'pasos', 5, 5);
registrar(cfg, '==== Etapa 15 lista: %s ====', dirR);
end

%% ======================================================================
function D = datos(conf, cfg, dirC, nw)
% Ventanas contiguas (etapa 13) con las familias pedidas + baselines + horas.
arch = fullfile(dirC, 'datos.mat');
if isfile(arch), S = load(arch); D = S.D; return, end
Cand = candidatos_rango(512, cfg);
baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for b = {'Re64a2', 'Re128a4'}
    c = Cand(strcmp({Cand.nombre}, b{1}));
    baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<AGROW>
end
D = struct();
for L = [512 1000]
    pre = sprintf('W%d', L);
    fams = unique(strsplit(strjoin(conf(startsWith(conf(:, 2), pre), 2), '+'), '+'));
    if isempty(fams{1}), continue, end
    [K, lay] = caract_registro(cfg, L, baseH);
    [~, ik] = ismember(fams, {K.nombre});
    S = load(fullfile(cfg.rutas.resultados, cfg.revision.dir_w, ['w_' pre '.mat']));
    [X, nom] = caract_desde_curvas(double(S.W.feat), K(ik), lay, cfg);
    T = removevars(S.W, {'feat', 'valida'}); T.feat = X;
    D.(pre) = T; D.([pre '_nombres']) = nom;
end
% baselines sobre la misma rejilla de 4 s (inicio cada 1000 muestras desde 0.5 s)
T = D.W1000; L = 1000; fs = cfg.fs;
[b8, a8] = butter(4, [8 30] / (fs/2));
ses = unique(T.sesion, 'stable');
partes = cell(numel(ses), 1);
dirS = fullfile(cfg.rutas.resultados, 'senales');
parfor (s = 1:numel(ses), nw)
    S = load(fullfile(dirS, [ses{s} '.mat']));
    filasS = find(strcmp(T.sesion, ses{s}));
    covs = zeros(16, 16, numel(filasS)); bp = zeros(numel(filasS), 48); hr = cell(numel(filasS), 1); q = 0;
    for j = 1:numel(S.X)
        [~, a] = fileparts(S.archivos{j});
        nf = nnz(strcmp(T.archivo(filasS), a));
        X = double(S.X{j}); Xf = filtfilt(b8, a8, X.').';
        ini = round(0.5 * fs) + 1 : L : size(X, 2) - L + 1;
        assert(numel(ini) == nf, 'etapa15: ventanas no alineadas en %s', a);
        for i = ini
            q = q + 1;
            V = Xf(:, i:i+L-1); covs(:, :, q) = V * V.' / L;
            [P, f] = pwelch(X(:, i:i+L-1).', hann(256), 128, 256, fs);
            bp(q, :) = [log(sum(P(f >= 4 & f < 8, :))), log(sum(P(f >= 8 & f < 13, :))), log(sum(P(f >= 13 & f < 30, :)))];
            hr{q} = S.archivos{j};
        end
    end
    partes{s} = struct('covs', covs, 'bp', bp, 'ruta', {hr});
end
covs = cellfun(@(p) p.covs, partes, 'UniformOutput', false);
D.covs = cat(3, covs{:});
D.bp = cell2mat(cellfun(@(p) p.bp, partes, 'UniformOutput', false));
rutas = vertcat(partes{:}); rutas = vertcat(rutas.ruta);
ch = cfg.canales.nombres;
D.BANDAS_nombres = [strcat('THETA_', ch), strcat('MU_', ch), strcat('BETA_', ch)];
% hora real de cada registro (fin de la grabación)
[ur, ~, ir] = unique(rutas);
h = NaT(numel(ur), 1);
for i = 1:numel(ur), h(i) = hora_registro(ur{i}); end
[~, arch] = cellfun(@fileparts, ur, 'UniformOutput', false);
D.horas = table(arch, h, 'VariableNames', {'archivo','hora'});
D.fila_archivo = arch(ir);
save(fullfile(dirC, 'datos.mat'), 'D', '-v7.3');
end

function c = campo(conj)
if startsWith(conj, 'W512'), c = 'W512'; else, c = 'W1000'; end
end

function [nom, ff] = config_feat(conj, D, T, cfg)
% Nombres de atributos y, para los baselines, función de características por fold.
ff = [];
switch conj
    case 'BANDAS'
        nom = D.BANDAS_nombres;
        bp = D.bp(fila_en_D(T, D), :);
        ff = @(tr, y) bp;
    case 'CSP'
        nom = arrayfun(@(i) sprintf('CSP%d', i), 1:6, 'UniformOutput', false);
        C = D.covs(:, :, fila_en_D(T, D));
        ff = @(tr, y) csp_feat(C, tr, y);
    otherwise
        Dn = D.([campo(conj) '_nombres']);
        pref = cellfun(@(n) n(1:find(n == '_', 1, 'last') - 1), Dn, 'UniformOutput', false);
        cols = ismember(pref, strsplit(conj, '+'));
        nom = Dn(cols);
        T0 = D.(campo(conj));
        X = T0.feat(fila_en(T, T0), cols);
        ff = @(tr, y) X;
end
end

function i = fila_en_D(T, D)
i = fila_en(T, D.W1000);
end

function i = fila_en(T, T0)
% Índice de cada fila de T en T0 (misma sesión, archivo y minuto).
k0 = strcat(T0.sesion, '|', T0.archivo, '|', cellstr(string(round(T0.t_archivo_min * 1e4))));
k = strcat(T.sesion, '|', T.archivo, '|', cellstr(string(round(T.t_archivo_min * 1e4))));
[ok, i] = ismember(k, k0);
assert(all(ok), 'etapa15: filas no encontradas');
end

function X = csp_feat(C, tr, y)
% CSP con las ventanas de entrenamiento: 3 pares de filtros, log-varianza normalizada.
n = size(C, 3);
Cn = C ./ reshape(arrayfun(@(i) trace(C(:, :, i)), 1:n), 1, 1, []);
Ca = mean(Cn(:, :, tr & y == 10), 3); Cb = mean(Cn(:, :, tr & y == 40), 3);
[Wf, Dd] = eig(Ca, Ca + Cb);
[~, o] = sort(diag(Dd)); Wf = Wf(:, o([1:3, end-2:end]));
X = zeros(n, 6);
for i = 1:n
    v = diag(Wf.' * C(:, :, i) * Wf);
    X(i, :) = log(v / sum(v));
end
end

function [a, ses] = acc_sesion(T, pred)
ses = unique(T.sesion, 'stable');
a = nan(numel(ses), 1);
for s = 1:numel(ses)
    k = strcmp(T.sesion, ses{s}) & isfinite(pred);
    if numel(unique(T.nivel(k))) == 2, a(s) = 100 * mean(pred(k) == T.nivel(k)); end
end
ok = isfinite(a); a = a(ok); ses = ses(ok);
end

function a = acc_sesion_vol(T, pred, vol)
y = T.nivel; v = ismember(T.sesion, vol); y(v) = 50 - y(v);
T.nivel = y;
a = acc_sesion(T, pred);
end

function ci = bootstrap_pac(a, pac)
% IC 95 % de la media por sesión remuestreando PACIENTES (bootstrap por grupos).
rs = RandStream('mt19937ar', 'Seed', 7);
up = unique(pac);
B = 10000; m = zeros(B, 1);
for b = 1:B
    sel = up(randi(rs, numel(up), numel(up), 1));
    v = cell2mat(cellfun(@(p) a(strcmp(pac, p)), sel, 'UniformOutput', false));
    m(b) = mean(v);
end
ci = quantile(m, [0.025 0.975]);
end

function [p, n] = perm_secuencial(T, nom, m, cfg, ff, obs, nw, archAv)
% Permutación de etiquetas por sesión completa (LOSO), con paro de Besag-Clifford.
sesC = unique(T.sesion);
nmax = 1000; h = 20; nula = nan(nmax, 1); n = 0;
while n < nmax
    idx = n+1 : min(nmax, n + 2 * max(nw, 1));
    parte = zeros(numel(idx), 1);
    parfor (j = 1:numel(idx), nw)
        rs = RandStream('mt19937ar', 'Seed', 500 + idx(j));
        vol = sesC(rand(rs, numel(sesC), 1) < 0.5);
        pq = predicciones_grupo(T, nom, m, cfg, T.sesion, vol, [], ff);
        parte(j) = mean(acc_sesion_vol(T, pq, vol));
    end
    nula(idx) = parte; n = idx(end);
    fid = fopen(archAv, 'w'); if fid > 0, fprintf(fid, '%d %d', n, nmax); fclose(fid); end
    if nnz(nula(1:n) >= obs) >= h, break, end
end
e = nnz(nula(1:n) >= obs);
if e >= h, p = e / n; else, p = (1 + e) / (1 + n); end
end

function [A, S, V] = votar(T, pred, nombre, esquema)
% Voto mayoritario por archivo y decisión por sesión (¿cuál registro es el de 40?).
[g, ses, arch] = findgroups(T.sesion, T.archivo);
f40 = splitapply(@(p) mean(p(isfinite(p)) == 40), pred, g);
niv = splitapply(@(x) x(1), T.nivel, g);
ok = isfinite(f40);
acierto = (f40 > 0.5) == (niv == 40);
A.acc = 100 * mean(acierto(ok)); A.p = 1 - binocdf(nnz(acierto(ok)) - 1, nnz(ok), 0.5);
us = unique(ses, 'stable'); dec = nan(numel(us), 1);
for s = 1:numel(us)
    i10 = strcmp(ses, us{s}) & niv == 10; i40 = strcmp(ses, us{s}) & niv == 40;
    if any(i10) && any(i40) && all(ok(i10 | i40)), dec(s) = f40(i40) > f40(i10); end
end
d = dec(isfinite(dec));
S.acc = 100 * mean(d); S.p = 1 - binocdf(nnz(d) - 1, numel(d), 0.5);
V = table(repmat({nombre}, numel(ses), 1), repmat({esquema}, numel(ses), 1), ses, arch, niv, f40, ...
    'VariableNames', {'configuracion','esquema','sesion','archivo','nivel','frac_pred40'});
end

function c = dias_completos(T, dia)
c = {};
for d = unique(dia).'
    k = strcmp(dia, d{1});
    if all(ismember({'KMI','FIS'}, unique(T.condicion(k)))) && numel(unique(T.archivo(k))) == 4, c{end+1} = d{1}; end %#ok<AGROW>
end
end

function O = fraccion_por_registro(T, pr, dia, horas, nombre)
[g, d, a, cnd] = findgroups(dia, T.archivo, T.condicion);
f40 = splitapply(@(p) mean(p == 40), pr, g);
niv = splitapply(@(x) x(1), T.nivel, g);
[~, ih] = ismember(a, horas.archivo);
h = horas.hora(ih);
minuto = zeros(numel(d), 1); pos = zeros(numel(d), 1);
for x = unique(d).'
    k = strcmp(d, x{1});
    minuto(k) = minutes(h(k) - min(h(k)));
    [~, o] = sort(h(k)); r = zeros(nnz(k), 1); r(o) = 1:nnz(k); pos(k) = r;
end
O = table(repmat({nombre}, numel(d), 1), d, a, cnd, niv, minuto, pos, f40, ...
    'VariableNames', {'configuracion','dia','archivo','condicion','nivel','minuto','posicion','frac_pred40'});
end

function [fila, txt] = modelo_mixto(O, nombre)
O.nivel40 = double(O.nivel == 40); O.fis = double(strcmp(O.condicion, 'FIS')); O.minuto10 = O.minuto / 10;
lme = fitlme(O, 'frac_pred40 ~ nivel40 + minuto10 + fis + (1|dia)');
c = lme.Coefficients;
b = @(n) c.Estimate(strcmp(c.Name, n)); pv = @(n) c.pValue(strcmp(c.Name, n));
con = @(c1, n1, c2, n2) mean(dif(O, c1, n1, c2, n2), 'omitnan');
fila = {nombre, b('nivel40'), pv('nivel40'), b('minuto10'), pv('minuto10'), b('fis'), pv('fis'), ...
    con('KMI', 10, 'FIS', 10), con('KMI', 40, 'FIS', 40), con('FIS', 40, 'FIS', 10), con('KMI', 40, 'KMI', 10), numel(unique(O.dia))};
txt = sprintf('nivel40 b = %+.3f (p = %.4f) | 10 min b = %+.3f (p = %.4f) | FIS b = %+.3f (p = %.3f) | KMI10-FIS10 %+.3f, FIS40-FIS10 %+.3f, KMI40-KMI10 %+.3f', ...
    fila{2:7}, fila{8}, fila{10}, fila{11});
end

function v = dif(O, c1, n1, c2, n2)
dias = unique(O.dia); v = nan(numel(dias), 1);
for i = 1:numel(dias)
    a = O.frac_pred40(strcmp(O.dia, dias{i}) & strcmp(O.condicion, c1) & O.nivel == n1);
    b = O.frac_pred40(strcmp(O.dia, dias{i}) & strcmp(O.condicion, c2) & O.nivel == n2);
    if ~isempty(a) && ~isempty(b), v(i) = a(1) - b(1); end
end
end

function r = ifelse(c, a, b)
if c, r = a; else, r = b; end
end

function avance(dirR, nombre, hecho, total)
fid = fopen(fullfile(dirR, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end
