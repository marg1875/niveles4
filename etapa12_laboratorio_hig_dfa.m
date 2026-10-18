function etapa12_laboratorio_hig_dfa(cfg)
%ETAPA12_LABORATORIO_HIG_DFA  Exploración amplia de Higuchi y DFA ("laboratorio").
%   Extiende la etapa 10 en tres ejes que no se habían probado:
%     1) Longitud de ventana: cfg.lab.L (p. ej. 1, 1.5, 2, 2.8 y 4 s).
%     2) DFA: orden de destendencia 1/2/3 y ventanas internas con y sin
%        50 % de traslape; escalas hasta L/3 (L/2 con traslape).
%     3) Higuchi: rango de k con k_max <= L/4.
%   Ventanas:
%     3c intra  mismas anclas de la etapa 3 (3 por bloque, cada 512 muestras).
%               Para L > 512 solo las anclas 1 y 3, para que las ventanas
%               vecinas no se traslapen (si no, la CV queda inflada).
%     contiguas para la LOSO y el excedente: inicio cada max(L, 512) muestras
%               (mismo número de ventanas para L <= 512; comparación justa).
%   Criterios por candidato (todos guardan el resultado por sesión):
%     A  3c     KMI 3 clases intra-sesión (13 sesiones x kNN/SVM/MLP ajustados)
%     B  exc    excedente: CV intra 10 vs 40 menos el control de mitades del
%               mismo archivo (descuenta la huella / deriva temporal), SVM ig8
%     C  LOSO   KMI 10 vs 40 dejando una sesión fuera, SVM ig8
%   Validación honesta: "anidado estricto" de la calibración completa. Para
%   cada sesión s se elige el candidato con el criterio A o B calculado SIN s
%   y se toma su exactitud LOSO en s (entrenada sin s). Como C ya guarda la
%   exactitud por sesión, no hay que reentrenar.
%   Combinaciones: finalistas por familia, L y criterio, pares HIG + DFA
%   dentro de la misma L y entre L distintas (multirresolución, ventanas
%   centradas en los mismos instantes), con los mismos criterios y anidado.
%   Final: evaluar_generalizacion y permutación secuencial de los ganadores
%   de las reglas anidadas.
%   -> resultados/laboratorio/*.csv, figuras/

lab = cfg.lab;
dirL = fullfile(cfg.rutas.resultados, lab.carpeta);
dirC = fullfile(dirL, 'cache');
dirFg = fullfile(dirL, 'figuras');
for d = {dirL, dirC, dirFg}
    if ~isfolder(d{1}), mkdir(d{1}); end
end
nw = usar_paralelo(cfg);
ch = cfg.canales.nombres; C = numel(ch);
Cand = candidatos_rango(512, cfg);
baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
for b = cfg.calib_hd.base_hurst
    c = Cand(strcmp({Cand.nombre}, b{1}));
    baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<AGROW>
end
cl = cfg; cl.calib_hd.hig_kmin = lab.hig_kmin; cl.calib_hd.hig_kmax = lab.hig_kmax;
cl.calib_hd.dfa_ns = lab.dfa_ns; cl.calib_hd.dfa_ordenes = lab.dfa_ordenes; cl.calib_hd.dfa_pasos = lab.dfa_pasos;
registrar(cfg, '==== Etapa 12: laboratorio Higuchi/DFA, L = %s, %d workers ====', mat2str(lab.L), nw);
nL = numel(lab.L);
avance(dirL, 'pasos', 0, 4 + nL);

%% Modelos
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste', 'hiperparametros_elegidos.csv'), 'TextType', 'char', 'Delimiter', ',');
H.Properties.VariableNames = {'modelo', 'especificacion'};
H = H(~strcmp(H.modelo, 'RF (ajustado)'), :);
mods3c = arrayfun(@(i) modelo_desde_id(H.especificacion{i}, 'ninguna'), 1:height(H));
etiq3c = strtok(H.modelo).';
Hc = readtable(fullfile(cfg.rutas.raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
mSVM = modelo_desde_id(Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')});

%% Sesiones 3c (metadatos de la etapa 3) y sesiones completas
dirArff = fullfile(cfg.rutas.resultados, 'arff', 'd512');
I = readtable(fullfile(dirArff, 'indice.csv'), 'TextType', 'char');
Ii = I(strcmp(I.condicion, 'KMI') & strcmp(I.esquema, '3c') & strcmp(I.alcance, 'intra'), :);
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
ses = ses(logical(ses.completa), :);
dirS = fullfile(cfg.rutas.resultados, 'senales');

%% 1) Por cada L: curvas, candidatos y los tres criterios
todos = table();
for iL = 1:nL
    L = lab.L(iL);
    pre = sprintf('W%d', L);
    cL = cl; cL.calib_hd.prefijo = pre;
    [K, lay] = candidatos_hig_dfa(cL, baseH, L);
    extractor = @(V) curvas_ventana(V, cL, baseH, lay);
    % a) curvas 3c intra
    arch3 = fullfile(dirC, sprintf('c3_%s.mat', pre));
    if ~isfile(arch3)
        t0 = tic; D3 = cell(height(Ii), 1);
        parfor (s = 1:height(Ii), nw)
            meta = readtable(fullfile(dirArff, strrep(Ii.archivo{s}, '.arff', '_meta.csv')), 'TextType', 'char');
            if L > 512, meta = meta(ismember(meta.ancla, [1 3]), :); end
            D3{s} = struct('id', Ii.id{s}, 'M', curvas_anclas(fullfile(dirS, [Ii.id{s} '.mat']), meta, extractor, L, cfg), ...
                'y', {arrayfun(@(c) sprintf('N%02d', c), meta.clase, 'UniformOutput', false)}, 'meta', meta);
        end
        save(arch3, 'D3', '-v7.3');
        registrar(cfg, '  %s: curvas 3c (%d ventanas) en %.1f min', pre, sum(cellfun(@(d) numel(d.y), D3)), toc(t0)/60);
    end
    % b) curvas en ventanas contiguas (KMI y FIS)
    archW = fullfile(dirC, sprintf('w_%s.mat', pre));
    if ~isfile(archW)
        t0 = tic; partes = cell(height(ses), 1); sid = ses.sesion; cnd = ses.condicion;
        rech = struct('id', 'R0', 'umbral_z', Inf, 'min_canales', 1);
        parfor (s = 1:height(ses), nw)
            partes{s} = ventanas_tiempo_sesion(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, baseH, L, cfg, ...
                rech, extractor, max(L, 512));
        end
        W = vertcat(partes{:});
        W.feat = matriz_caracteristicas(W);
        save(archW, 'W', '-v7.3');
        registrar(cfg, '  %s: %d ventanas contiguas en %.1f min', pre, height(W), toc(t0)/60);
    end
    % c) criterio A (3c intra), por sesión
    archA = fullfile(dirC, sprintf('A_%s.mat', pre));
    if ~isfile(archA)
        t0 = tic;
        S = load(arch3); D3 = S.D3;
        RA = cell(numel(D3), 1);
        parfor (s = 1:numel(D3), nw)
            RA{s} = criterio_3c(D3{s}, K, lay, cL, mods3c, etiq3c);
        end
        A = vertcat(RA{:});
        save(archA, 'A');
        registrar(cfg, '  %s: criterio 3c, %d candidatos (%.1f min)', pre, numel(K), toc(t0)/60);
    end
    % d) criterios B (excedente) y C (LOSO), por bloques de candidatos
    archBC = fullfile(dirC, sprintf('BC_%s.mat', pre));
    if ~isfile(archBC)
        t0 = tic;
        S = load(archW); WK = S.W(strcmp(S.W.condicion, 'KMI'), :);
        Mk = parallel.pool.Constant(WK.feat);
        meta = removevars(WK, 'feat');
        bloques = num2cell(reshape([1:numel(K), nan(1, mod(-numel(K), lab.bloque))], lab.bloque, []), 1);
        bloques = cellfun(@(b) b(isfinite(b)), bloques, 'UniformOutput', false);
        archB = cellfun(@(b) fullfile(dirC, sprintf('bc_%s_%04d.mat', pre, b(1))), bloques, 'UniformOutput', false);
        pend = find(~cellfun(@isfile, archB));
        hechos = numel(bloques) - numel(pend);
        for g0 = 1:nw:numel(pend)
            grupo = pend(g0:min(numel(pend), g0 + nw - 1));
            parfor (q = 1:numel(grupo), nw)
                b = grupo(q);
                R = criterios_bc(Mk.Value, meta, K(bloques{b}), lay, cL, mSVM);
                guardar(archB{b}, R);
            end
            hechos = hechos + numel(grupo);
            avance(dirL, ['bc_' pre], hechos, numel(bloques));
        end
        BC = vertcat(cellfun(@(a) getfield(load(a), 'R'), archB, 'UniformOutput', false)); %#ok<GFLD>
        BC = vertcat(BC{:});
        save(archBC, 'BC');
        registrar(cfg, '  %s: criterios excedente y LOSO, %d candidatos (%.1f min)', pre, numel(K), toc(t0)/60);
    end
    S = load(archA); A = S.A; S = load(archBC); BC = S.BC;
    todos = [todos; resumen_L(A, BC, K, L)]; %#ok<AGROW>
    avance(dirL, 'pasos', iL, 4 + nL);
end
writetable(todos, fullfile(dirL, 'candidatos_resumen.csv'));

%% 2) Rankings y anidado estricto (candidatos individuales)
[A, BC] = cargar_todos(dirC, lab.L);
reglas = anidado_estricto(A, BC, todos, {'todos', 'HIG', 'DFA'});
for iL = 1:nL
    reglas = [reglas; anidado_estricto(A, BC, todos(todos.L == lab.L(iL), :), {sprintf('W%d', lab.L(iL))})]; %#ok<AGROW>
end
writetable(reglas, fullfile(dirL, 'anidado_individual.csv'));
for i = 1:height(reglas)
    registrar(cfg, '  anidado %-6s regla %-4s: LOSO %.1f %% (%d/%d > 50, p_W = %.3f) | más elegido: %s', ...
        reglas.pool{i}, reglas.regla{i}, reglas.acc(i), reglas.sobre50(i), reglas.n(i), reglas.p_wilcoxon(i), reglas.mas_elegido{i});
end
graficar(todos, lab.L, dirFg);
avance(dirL, 'pasos', nL + 1, 4 + nL);

%% 3) Combinaciones HIG + DFA (misma L y multirresolución)
fin = finalistas(todos, lab.top);
pares = {};
for a = 1:numel(fin.HIG), for b = 1:numel(fin.DFA), pares{end+1} = [fin.HIG{a} '+' fin.DFA{b}]; end, end %#ok<AGROW>
pares = [pares, refs_combo(todos)];
pares = unique(pares, 'stable');
registrar(cfg, '  combinaciones: %d pares (finalistas HIG %d x DFA %d)', numel(pares), numel(fin.HIG), numel(fin.DFA));
Pc = combinaciones(pares, lab, dirC, dirL, baseH, cl, cfg, mods3c, etiq3c, mSVM, Ii, ses, dirS, dirArff, nw);
writetable(Pc.resumen, fullfile(dirL, 'combinaciones_resumen.csv'));
reglasC = anidado_estricto(Pc.A, Pc.BC, Pc.resumen, {'combos'});
reglasC = [reglasC; anidado_estricto([A; Pc.A], [BC; Pc.BC], [todos; Pc.resumen], {'todo'})];
writetable(reglasC, fullfile(dirL, 'anidado_combinaciones.csv'));
for i = 1:height(reglasC)
    registrar(cfg, '  anidado %-6s regla %-4s: LOSO %.1f %% (%d/%d > 50, p_W = %.3f) | más elegido: %s', ...
        reglasC.pool{i}, reglasC.regla{i}, reglasC.acc(i), reglasC.sobre50(i), reglasC.n(i), reglasC.p_wilcoxon(i), reglasC.mas_elegido{i});
end
avance(dirL, 'pasos', nL + 2, 4 + nL);

%% 4) Validación final de los ganadores: generalización completa + permutación secuencial
R = [todos; Pc.resumen];
gan = unique([cellfun(@(c) c, ganadores(R), 'UniformOutput', false), {'W512HIG3a8+W512DFAc4a24', 'W512Re64a2+W512Re128a4'}], 'stable');
gan = gan(ismember(gan, R.conjunto));
F = validacion_final(gan, lab, dirC, dirL, baseH, cl, cfg, mSVM, ses, dirS, nw);
writetable(F, fullfile(dirL, 'validacion_final.csv'));
avance(dirL, 'pasos', nL + 4, 4 + nL);
registrar(cfg, '==== Etapa 12 lista: %s ====', dirL);
end

%% ======================================================================
function M = curvas_anclas(rutaSenal, meta, extractor, L, cfg)
% Curvas de las ventanas centradas en las anclas de la etapa 3 (meta).
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

function R = criterio_3c(D, K, lay, cL, mods, etiq)
% Criterio A: exactitud CV 5-fold de cada candidato en una sesión 3c.
weka_iniciar(cL);
X = caracteristicas_desde_curvas(D.M, K, lay, cL);
C = lay.C;
acc = nan(numel(K), numel(mods));
for k = 1:numel(K)
    cols = (k-1)*C + (1:C);
    if all(isnan(X(:, cols)), 'all'), continue, end
    Dk = weka_instancias(D.id, arrayfun(@(i) sprintf('a%d', i), 1:C, 'UniformOutput', false), X(:, cols), {'N00','N10','N40'}, D.y);
    for m = 1:numel(mods)
        acc(k, m) = weka_validar(Dk, weka_clasificador_param(mods(m), cL), cL).acc;
    end
end
R = table(repmat({D.id}, numel(K), 1), {K.nombre}.', mean(acc, 2), acc, ...
    'VariableNames', {'sesion3c','conjunto','acc3c','acc3c_clf'});
end

function R = criterios_bc(M, meta, K, lay, cL, m)
% Criterios B (excedente intra 10/40 sobre mitades) y C (LOSO) por sesión.
weka_iniciar(cL);
X = caracteristicas_desde_curvas(M, K, lay, cL);
C = lay.C;
ses = unique(meta.sesion, 'stable');
filas = cell(numel(K), 1);
for k = 1:numel(K)
    cols = (k-1)*C + (1:C);
    nom = strcat(K(k).nombre, '_', cL.canales.nombres);
    T = meta; T.feat = X(:, cols);
    if all(isnan(T.feat), 'all'), filas{k} = table(); continue, end
    loso = loso_niveles_todas(T, nom, m, cL);
    exc = nan(numel(ses), 1); intra = exc; mitad = exc;
    for s = 1:numel(ses)
        Ts = T(strcmp(T.sesion, ses{s}), :);
        intra(s) = cv_mem(Ts, nom, m, cL, 'nivel');
        arch = unique(Ts.archivo);
        mitad(s) = mean(cellfun(@(a) cv_mem(Ts(strcmp(Ts.archivo, a), :), nom, m, cL, 'mitad'), arch), 'omitnan');
        exc(s) = intra(s) - mitad(s);
    end
    filas{k} = table(repmat({K(k).nombre}, numel(ses), 1), ses, loso, intra, mitad, exc, ...
        'VariableNames', {'conjunto','sesion','loso','intra','mitad','exc'});
end
R = vertcat(filas{:});
end

function acc = loso_niveles_todas(T, nom, m, cfg)
% LOSO con una fila por sesión de T (NaN si la sesión no es evaluable).
ses = unique(T.sesion, 'stable');
acc = nan(numel(ses), 1);
X = matriz_caracteristicas(T);
for s = 1:numel(ses)
    i = strcmp(T.sesion, ses{s});
    X(i, :) = (X(i, :) - mean(X(i, :), 1, 'omitnan')) ./ std(X(i, :), 0, 1, 'omitnan');
end
y = arrayfun(@(n) sprintf('N%02d', n), T.nivel, 'UniformOutput', false);
for s = 1:numel(ses)
    te = strcmp(T.sesion, ses{s});
    if numel(unique(y(te))) < 2, continue, end
    Dtr = weka_instancias('tr', nom, X(~te, :), {'N10','N40'}, y(~te));
    Dte = weka_instancias('te', nom, X(te, :), {'N10','N40'}, y(te));
    cls = weka_clasificador_param(m, cfg); cls.buildClassifier(Dtr);
    ok = 0;
    for i = 0:Dte.numInstances()-1
        ok = ok + (cls.classifyInstance(Dte.instance(i)) == Dte.instance(i).classValue());
    end
    acc(s) = 100 * ok / Dte.numInstances();
end
end

function a = cv_mem(T, nom, m, cfg, tipo)
% CV 5-fold en memoria dentro de una sesión (nivel) o de un archivo (mitades).
if strcmp(tipo, 'nivel')
    y = arrayfun(@(n) sprintf('N%02d', n), T.nivel, 'UniformOutput', false); cl = {'N10','N40'};
else
    y = repmat({'M1'}, height(T), 1); y(T.t_archivo_min > median(T.t_archivo_min)) = {'M2'}; cl = {'M1','M2'};
end
a = NaN;
if min(cellfun(@(c) nnz(strcmp(y, c)), cl)) < cfg.weka.folds, return, end
D = weka_instancias('cv', nom, matriz_caracteristicas(T), cl, y);
a = weka_validar(D, weka_clasificador_param(m, cfg), cfg).acc;
end

function Rs = resumen_L(A, BC, K, L)
% Una fila por candidato: medias de los tres criterios.
[g, cj] = findgroups(A.conjunto);
a3 = splitapply(@(x) mean(x, 'omitnan'), A.acc3c, g);
[g2, cj2] = findgroups(BC.conjunto);
lo = splitapply(@(x) mean(x, 'omitnan'), BC.loso, g2);
ex = splitapply(@(x) mean(x, 'omitnan'), BC.exc, g2);
s50 = splitapply(@(x) nnz(x > 50), BC.loso, g2);
Rs = table(cj, repmat(L, numel(cj), 1), 'VariableNames', {'conjunto','L'});
[~, ia] = ismember(Rs.conjunto, {K.nombre});
tip = repmat({'COMB'}, height(Rs), 1); tip(ia > 0) = {K(ia(ia > 0)).tipo};
Rs.familia = tip;
Rs.acc3c = a3;
[~, ib] = ismember(Rs.conjunto, cj2);
Rs.loso = nan(height(Rs), 1); Rs.exc = Rs.loso; Rs.sobre50 = Rs.loso;
Rs.loso(ib > 0) = lo(ib(ib > 0)); Rs.exc(ib > 0) = ex(ib(ib > 0)); Rs.sobre50(ib > 0) = s50(ib(ib > 0));
end

function [A, BC] = cargar_todos(dirC, Ls)
A = table(); BC = table();
for L = Ls
    S = load(fullfile(dirC, sprintf('A_W%d.mat', L))); A = [A; S.A]; %#ok<AGROW>
    S = load(fullfile(dirC, sprintf('BC_W%d.mat', L))); BC = [BC; S.BC]; %#ok<AGROW>
end
end

function T = anidado_estricto(A, BC, resumen, pools)
% Para cada sesión de prueba s, elige el candidato con el criterio (3c o
% excedente) calculado sin s y toma su LOSO en s. También reporta la regla
% "LOSO" (elegir por la media de LOSO sin s), que es solo aproximadamente
% honesta: los modelos de las demás sesiones sí vieron s al entrenar.
filas = {};
sesL = unique(BC.sesion, 'stable');
A.sesion = cellfun(@(x) mapear(x), A.sesion3c, 'UniformOutput', false);
for p = pools
    switch p{1}
        case {'todos', 'todo', 'combos'}, cand = resumen.conjunto;
        case {'HIG', 'DFA'}, cand = resumen.conjunto(strcmp(resumen.familia, p{1}));
        otherwise, cand = resumen.conjunto(startsWith(resumen.conjunto, p{1}));
    end
    Ac = A(ismember(A.conjunto, cand), :);
    Bc = BC(ismember(BC.conjunto, cand), :);
    for regla = {'3c', 'exc', 'LOSO'}
        acc = nan(numel(sesL), 1); eleg = cell(numel(sesL), 1);
        for s = 1:numel(sesL)
            switch regla{1}
                case '3c',  T0 = Ac(~strcmp(Ac.sesion, sesL{s}), :); v = T0.acc3c;
                case 'exc', T0 = Bc(~strcmp(Bc.sesion, sesL{s}), :); v = T0.exc;
                case 'LOSO', T0 = Bc(~strcmp(Bc.sesion, sesL{s}), :); v = T0.loso;
            end
            [g, cj] = findgroups(T0.conjunto);
            [~, j] = max(splitapply(@(x) mean(x, 'omitnan'), v, g));
            eleg{s} = cj{j};
            r = Bc.loso(strcmp(Bc.conjunto, cj{j}) & strcmp(Bc.sesion, sesL{s}));
            if ~isempty(r), acc(s) = r(1); end
        end
        ok = isfinite(acc);
        W = wilcoxon_holm([acc(ok), 50 * ones(nnz(ok), 1)], {'a', 'b'});
        [u, ~, ic] = unique(eleg(ok)); [~, im] = max(accumarray(ic, 1));
        filas(end+1, :) = {p{1}, regla{1}, mean(acc(ok)), nnz(acc(ok) > 50), nnz(ok), W.p, u{im}, strjoin(eleg(ok).', ' | ')}; %#ok<AGROW>
    end
end
T = cell2table(filas, 'VariableNames', {'pool','regla','acc','sobre50','n','p_wilcoxon','mas_elegido','elegidos'});
end

function s = mapear(id3c)
% 'KMI_3c_Px006_M1' -> 'Px006_M1_KMI'
t = regexp(id3c, '(Px\d+)_M(\d)', 'tokens', 'once');
s = sprintf('%s_M%s_KMI', t{1}, t{2});
end

function f = finalistas(R, n)
% Mejores n por familia según 3c, excedente y LOSO (unión, todas las L).
f = struct('HIG', {{}}, 'DFA', {{}});
for fam = {'HIG', 'DFA'}
    Rf = R(strcmp(R.familia, fam{1}), :);
    c = {};
    for crit = {'acc3c', 'exc', 'loso'}
        Rs = sortrows(Rf, crit{1}, 'descend');
        c = [c; Rs.conjunto(1:min(n, height(Rs)))]; %#ok<AGROW>
    end
    f.(fam{1}) = unique(c, 'stable').';
end
end

function c = refs_combo(R)
% Combinaciones de referencia que existen en la rejilla.
c = {};
for x = {'W512HIG3a8+W512DFAc4a24', 'W512Re64a2+W512Re128a4', 'W512Re16a2+W512Re64a4'}
    p = strsplit(x{1}, '+');
    if all(ismember(p, R.conjunto)), c{end+1} = x{1}; end %#ok<AGROW>
end
end

function P = combinaciones(pares, lab, dirC, dirL, baseH, cl, cfg, mods3c, etiq3c, mSVM, Ii, ses, dirS, dirArff, nw)
% Pares (posiblemente de L distintas): mismas anclas para 3c (intersección
% de anclas si alguna L > 512) y ventanas contiguas centradas en los mismos
% instantes (paso max(L) del par) para LOSO y excedente.
arch = fullfile(dirC, 'combinaciones.mat');
if isfile(arch), S = load(arch); P = S.P; return, end
t0 = tic;
Ls = unique(cellfun(@(c) str2double(regexp(c, '^W(\d+)', 'tokens', 'once')), strsplit(strjoin(pares, '+'), '+')));
famsPorL = containers.Map('KeyType', 'double', 'ValueType', 'any');
for L = Ls
    famsPorL(L) = unique(regexp(strjoin(pares, '+'), sprintf('W%d[A-Za-z0-9]+', L), 'match'));
end
% --- ventanas centradas comunes: centros cada max(Ls) muestras
paso = max(max(Ls), 512);
partes = cell(height(ses), 1); sid = ses.sesion; cnd = ses.condicion;
parfor (s = 1:height(ses), nw)
    partes{s} = ventanas_centradas(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, Ls, paso, famsPorL, baseH, cl, cfg);
end
Wc = vertcat(partes{:});
% --- 3c: mismas anclas (anclas 1 y 3 si hay L > 512)
D3 = cell(height(Ii), 1);
parfor (s = 1:height(Ii), nw)
    meta = readtable(fullfile(dirArff, strrep(Ii.archivo{s}, '.arff', '_meta.csv')), 'TextType', 'char');
    if any(Ls > 512), meta = meta(ismember(meta.ancla, [1 3]), :); end
    D3{s} = anclas_multi(fullfile(dirS, [Ii.id{s} '.mat']), meta, Ls, famsPorL, baseH, cl, cfg, Ii.id{s});
end
nombres = partes{1}.Properties.UserData;
WK = Wc(strcmp(Wc.condicion, 'KMI'), :);
RA = cell(numel(D3), 1); RB = cell(numel(pares), 1);
parfor (s = 1:numel(D3), nw)
    RA{s} = criterio_3c_cols(D3{s}, pares, nombres, cfg, mods3c);
end
Xk = parallel.pool.Constant(matriz_caracteristicas(WK));
meta = removevars(WK, 'feat');
parfor (k = 1:numel(pares), nw)
    RB{k} = criterio_bc_cols(Xk.Value, meta, pares{k}, nombres, cfg, mSVM);
end
P.A = vertcat(RA{:}); P.BC = vertcat(RB{:});
P.resumen = resumen_L(P.A, P.BC, struct('nombre', {}, 'tipo', {}), 0);
P.resumen.familia(:) = {'COMB'};
save(arch, 'P');
registrar(cfg, '  combinaciones evaluadas: %d (%.1f min)', numel(pares), toc(t0)/60);
avance(dirL, 'combinaciones', 1, 1);
end

function T = ventanas_centradas(rutaS, sid, cond, Ls, paso, famsPorL, baseH, cl, cfg)
% Ventanas de cada L centradas en los mismos instantes; solo las familias pedidas.
S = load(rutaS);
nivel = S.niveles;
stimsMax = cellfun(@max, S.stims);
if numel(S.X) == 2 && all(sort(stimsMax) == [10 40]), nivel = stimsMax; end
filas = {}; nombres = {};
for j = 1:numel(S.X)
    X = double(S.X{j}); N = size(X, 2);
    [~, arch] = fileparts(S.archivos{j});
    centros = round(0.5 * cfg.fs) + ceil(max(Ls) / 2) : paso : N - ceil(max(Ls) / 2);
    for c = centros
        f = []; nombres = {};
        for L = Ls
            [x, n] = feats_L(X(:, c - floor(L/2) + (1:L)), L, famsPorL(L), baseH, cl);
            f = [f, x]; nombres = [nombres, n]; %#ok<AGROW>
        end
        filas(end+1, :) = {sid, cond, arch, nivel(j), c / cfg.fs / 60, f}; %#ok<AGROW>
    end
end
T = cell2table(filas, 'VariableNames', {'sesion','condicion','archivo','nivel','t_archivo_min','feat'});
T.Properties.UserData = nombres;
end

function D = anclas_multi(rutaSenal, meta, Ls, famsPorL, baseH, cl, cfg, id)
S = load(rutaSenal);
nombresArch = cellfun(@(r) regexprep(r, '^.*[\\/]|\.mat$', ''), S.archivos, 'UniformOutput', false);
X = [];
for j = 1:numel(S.X)
    filas = find(strcmp(meta.archivo, nombresArch{j}));
    if isempty(filas), continue, end
    A = segmentar_por_estimulos(size(S.X{j}, 2), S.stims{j}, cfg);
    for f = filas.'
        k = find(A.bloque == meta.bloque(f) & A.ancla == meta.ancla(f), 1);
        x = [];
        for L = Ls
            x = [x, feats_L(double(extraer_ventanas(S.X{j}, A.centro(k), L)), L, famsPorL(L), baseH, cl)]; %#ok<AGROW>
        end
        if isempty(X), X = nan(height(meta), numel(x)); end
        X(f, :) = x;
    end
end
D = struct('id', id, 'X', X, 'y', {arrayfun(@(c) sprintf('N%02d', c), meta.clase, 'UniformOutput', false)});
end

function [x, n] = feats_L(V, L, fams, baseH, cl)
% Familias 'fams' (con prefijo W<L>) de una ventana de longitud L.
cL = cl; cL.calib_hd.prefijo = sprintf('W%d', L);
[K, lay] = candidatos_hig_dfa(cL, baseH, L);
K = K(ismember({K.nombre}, fams));
[x, n] = caracteristicas_desde_curvas(curvas_ventana(V, cL, baseH, lay), K, lay, cL);
end

function R = criterio_3c_cols(D, pares, nombres, cfg, mods)
weka_iniciar(cfg);
acc = nan(numel(pares), numel(mods));
for k = 1:numel(pares)
    cols = columnas(nombres, pares{k});
    Dk = weka_instancias(D.id, nombres(cols), D.X(:, cols), {'N00','N10','N40'}, D.y);
    for m = 1:numel(mods), acc(k, m) = weka_validar(Dk, weka_clasificador_param(mods(m), cfg), cfg).acc; end
end
R = table(repmat({D.id}, numel(pares), 1), pares(:), mean(acc, 2), acc, ...
    'VariableNames', {'sesion3c','conjunto','acc3c','acc3c_clf'});
end

function R = criterio_bc_cols(X, meta, par, nombres, cfg, m)
weka_iniciar(cfg);
cols = columnas(nombres, par);
T = meta; T.feat = X(:, cols);
ses = unique(meta.sesion, 'stable');
loso = loso_niveles_todas(T, nombres(cols), m, cfg);
intra = nan(numel(ses), 1); mitad = intra;
for s = 1:numel(ses)
    Ts = T(strcmp(T.sesion, ses{s}), :);
    intra(s) = cv_mem(Ts, nombres(cols), m, cfg, 'nivel');
    mitad(s) = mean(cellfun(@(a) cv_mem(Ts(strcmp(Ts.archivo, a), :), nombres(cols), m, cfg, 'mitad'), unique(Ts.archivo)), 'omitnan');
end
R = table(repmat({par}, numel(ses), 1), ses, loso, intra, mitad, intra - mitad, ...
    'VariableNames', {'conjunto','sesion','loso','intra','mitad','exc'});
end

function g = ganadores(R)
% Mejor por criterio entre individuales y combinaciones.
g = {};
for crit = {'acc3c', 'exc', 'loso'}
    Rs = sortrows(R(isfinite(R.(crit{1})), :), crit{1}, 'descend');
    g{end+1} = Rs.conjunto{1}; %#ok<AGROW>
end
end

function F = validacion_final(gan, lab, dirC, dirL, baseH, cl, cfg, m, ses, dirS, nw)
% Métricas completas (KMI y FIS) y permutación secuencial de los ganadores.
Ls = unique(cellfun(@(c) str2double(regexp(c, '^W(\d+)', 'tokens', 'once')), strsplit(strjoin(gan, '+'), '+')));
famsPorL = containers.Map('KeyType', 'double', 'ValueType', 'any');
for L = Ls, famsPorL(L) = unique(regexp(strjoin(gan, '+'), sprintf('W%d[A-Za-z0-9]+', L), 'match')); end
archW = fullfile(dirC, 'ventanas_final.mat');
if isfile(archW)
    S = load(archW); W = S.W; nombres = S.nombres;
else
    partes = cell(height(ses), 1); sid = ses.sesion; cnd = ses.condicion; paso = max(max(Ls), 512);
    parfor (s = 1:height(ses), nw)
        partes{s} = ventanas_centradas(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, Ls, paso, famsPorL, baseH, cl, cfg);
    end
    W = vertcat(partes{:}); nombres = partes{1}.Properties.UserData;
    W.feat = matriz_caracteristicas(W);
    save(archW, 'W', 'nombres', '-v7.3');
end
filas = cell(numel(gan), 1);
for i = 1:numel(gan)
    archF = fullfile(dirC, sprintf('final_%02d.mat', i));
    if isfile(archF), S = load(archF); filas{i} = S.fila; continue, end
    t0 = tic;
    cols = columnas(nombres, gan{i});
    Tk = removevars(W, 'feat'); Tk.feat = W.feat(:, cols);
    G = evaluar_generalizacion(Tk, nombres(cols), m, cfg);
    TK = Tk(strcmp(Tk.condicion, 'KMI'), :);
    [obs, nula, p] = permutacion_secuencial(TK, nombres(cols), m, cfg, nw, 1000, 20, ...
        fullfile(dirL, sprintf('perm_%02d.avance.txt', i)));
    fila = [table(gan(i), numel(cols), obs, mean(nula), numel(nula), p, ...
        'VariableNames', {'conjunto','n_atributos','loso_kmi_perm','nula_media','n_perm','p_perm'}), G];
    save(archF, 'fila');
    filas{i} = fila;
    registrar(cfg, '  final %s: LOSO %.1f %% (p perm = %.4f, %d perm.) | FIS %.1f %% | excedente %.1f (%.1f min)', ...
        gan{i}, obs, p, numel(nula), G.loso_fis, G.excedente_kmi, toc(t0)/60);
end
F = vertcat(filas{:});
end

function cols = columnas(nombres, conjunto)
pref = cellfun(@(n) n(1:find(n == '_', 1, 'last') - 1), nombres, 'UniformOutput', false);
cols = find(ismember(pref, strsplit(conjunto, '+')));
end

function graficar(R, Ls, dirFg)
% Mejor valor por L y familia/variante en cada criterio.
f = figure('Visible', 'off', 'Position', [100 100 1400 420]);
crit = {'acc3c', 'KMI 3c intra (%)'; 'exc', 'Excedente sobre mitades (pts)'; 'loso', 'LOSO KMI 10 vs 40 (%)'};
var = {'HIG', '^W\d+HIG'; 'DFA1', '^W\d+DFAl\d'; 'DFA1o', '^W\d+DFAlo'; 'DFA2', '^W\d+DFAc\d'; ...
       'DFA2o', '^W\d+DFAco'; 'DFA3', '^W\d+DFAt\d'; 'DFA3o', '^W\d+DFAto'};
for q = 1:3
    subplot(1, 3, q); hold on;
    for v = 1:size(var, 1)
        y = nan(size(Ls));
        for i = 1:numel(Ls)
            s = R.L == Ls(i) & ~cellfun(@isempty, regexp(R.conjunto, var{v, 2}, 'once'));
            if any(s), y(i) = max(R.(crit{q, 1})(s)); end
        end
        plot(Ls / 250, y, '-o', 'DisplayName', var{v, 1});
    end
    xlabel('Ventana (s)'); title(crit{q, 2}); grid on;
    if q == 1, legend('Location', 'best'); end
end
exportgraphics(f, fullfile(dirFg, 'laboratorio_por_ventana.png'), 'Resolution', 130);
close(f);
end

function guardar(ruta, R)
tmp = [ruta '.tmp.mat']; save(tmp, 'R'); movefile(tmp, ruta, 'f');
end

function avance(dirL, nombre, hecho, total)
fid = fopen(fullfile(dirL, [nombre '.avance.txt']), 'w');
if fid > 0, fprintf(fid, '%d %d', hecho, total); fclose(fid); end
end
