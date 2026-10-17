function etapa8_reposo_eventos(cfg)
%ETAPA8_REPOSO_EVENTOS  Reposo y tarea redefinidos; pruebas del sesgo de archivo.
%   Motivo: analisis_confundido_reposo.m mostró que el reposo se clasificaba
%   como el nivel de SU archivo (el modelo reconocía la grabación). Además,
%   la sesión imaginada (KMI) fue ASÍNCRONA: su stims no indica cuándo
%   imaginó el paciente. Aquí:
%     - FIS (síncrona): tarea = bloques de stims (eventos_fis).
%     - KMI (asíncrona): tarea = intentos detectados por ERD (detectar_erd).
%       El detector se valida en FIS, donde los bloques son conocidos.
%     - Reposo, dos definiciones (misma longitud de ventana por dataset):
%         A: primeros 5 s del registro (ventanas de cfg.eventos.L_A = 1 s)
%         B: [fin + 0.5, fin + 0.8] s tras cada estímulo / intento (0.3 s)
%     - Variante Δ: a cada característica se le resta la media del reposo
%       inicial de su propio archivo (elimina la "huella" de la grabación).
%   Pruebas (intra-sujeto y global, CV 5-fold x cfg.eventos.semillas):
%     3c          reposo / 10 % / 40 %  (+ índice de herencia del reposo)
%     2c          10 % vs 40 %
%     huella      reposo del archivo 10 % vs reposo del archivo 40 %  (ideal ~50 %)
%     tarea_rep   tarea vs reposo (ambos archivos; sin confundido de archivo)
%   -> resultados/eventos/*.csv, figuras; ARFF en resultados/arff/reposo_*/

cfg.eventos.sufijo = '';
if cfg.modo_prueba, cfg.eventos.sufijo = '_prueba'; end
dirE = fullfile(cfg.rutas.resultados, ['eventos' cfg.eventos.sufijo]);
dirEC = fullfile(dirE, 'cache'); dirEF = fullfile(dirE, 'figuras');
for d = {dirE, dirEC, dirEF}, if ~isfolder(d{1}), mkdir(d{1}); end, end
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
ses = ses(logical(ses.completa), :);
if cfg.modo_prueba, ses = ses(strcmp(ses.paciente, 'Px006') & ses.mes == 1, :); end
nw = usar_paralelo(cfg);
registrar(cfg, 'Etapa 8: reposo/eventos en %d sesiones completas.', height(ses));

%% 1) Eventos, ventanas y características
archF = fullfile(dirE, 'caracteristicas_eventos.mat');
if ~isfile(archF)
    t0 = tic;
    dirS = fullfile(cfg.rutas.resultados, 'senales');
    partes = cell(height(ses), 1); evs = partes; vals = partes;
    sid = ses.sesion; cond = ses.condicion; pac = ses.paciente; mes = ses.mes;
    parfor (s = 1:height(ses), nw)
        [partes{s}, evs{s}, vals{s}] = procesar_sesion(fullfile(dirS, [sid{s} '.mat']), ...
            sid{s}, pac{s}, mes(s), cond{s}, cfg);
    end
    F = vertcat(partes{:});
    EV = vertcat(evs{:});
    VA = vertcat(vals{:});
    save(archF, 'F', 'EV', 'VA', '-v7.3');
    writetable(EV, fullfile(dirE, 'eventos_por_archivo.csv'));
    writetable(VA, fullfile(dirE, 'validacion_erd_fis.csv'));
    registrar(cfg, '  eventos y características listos (%.1f min).', toc(t0)/60);
end
S = load(archF); F = S.F; VA = S.VA; EV = S.EV;
if ~isempty(VA)
    registrar(cfg, '  ERD en FIS: %.1f %% de los eventos caen en bloques de estímulo (azar = %.1f %%); %.1f %% de bloques detectados.', ...
        100*mean(VA.prop_eventos_en_tarea, 'omitnan'), 100*mean(VA.frac_tiempo_tarea), 100*mean(VA.prop_bloques_detectados));
end
k = EV(strcmp(EV.condicion, 'KMI'), :);
registrar(cfg, '  KMI: %.1f intentos ERD por archivo (rango %d-%d).', mean(k.n_eventos), min(k.n_eventos), max(k.n_eventos));

%% 2) Datasets ARFF
nombresF = {};
for k2 = 1:numel(cfg.eventos.familias)
    nombresF = [nombresF, strcat(cfg.eventos.familias(k2).nombre, '_', cfg.canales.nombres)]; %#ok<AGROW>
end
unidades = construir_datasets(F, nombresF, cfg);
registrar(cfg, '  %d datasets ARFF.', numel(unidades));

%% 3) Evaluación
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), 'Delimiter', ',', 'TextType', 'char');
H.Properties.VariableNames = {'modelo','especificacion'};
modelos = [];
for i = find(~strcmp(H.modelo, 'RF (ajustado)')).'
    modelos = [modelos, modelo_desde_id(H.especificacion{i})]; %#ok<AGROW>
end
modelos = [modelos, especificar_modelo('kNN', struct()), especificar_modelo('SVM', struct()), ...
           especificar_modelo('MLP', struct())];
etiq = [strtok(H.modelo(~strcmp(H.modelo, 'RF (ajustado)'))).' + " ajustado", "kNN artículo", "SVM artículo", "MLP artículo"];
etiq = cellstr(etiq);
t0 = tic;
parfor (u = 1:numel(unidades), nw)
    evaluar_unidad(unidades(u), modelos, etiq, cfg, fullfile(dirEC, [unidades(u).clave '.mat']));
end
registrar(cfg, '  evaluación lista (%.1f min).', toc(t0)/60);

%% 4) Resumen
partes = cell(numel(unidades), 1);
for u = 1:numel(unidades)
    a = fullfile(dirEC, [unidades(u).clave '.mat']);
    if isfile(a)
        Sx = load(a, 'R');
        Rx = Sx.R;
        if iscell(Rx.azar)                           % cachés con filas de error mal formadas
            Rx.azar = cellfun(@(x) sum(double(x)) * isnumeric(x) + NaN * ~isnumeric(x), Rx.azar);
        end
        if height(Rx) > 0, partes{u} = Rx; end       % los datasets sin filas se omiten
    end
end
R = vertcat(partes{:});
writetable(R, fullfile(dirE, 'resultados_por_dataset.csv'));
T = resumen_final(R);
writetable(T, fullfile(dirE, 'tabla_final_reposo.csv'));
graficar_eventos(T, dirEF);
registrar(cfg, 'Etapa 8 lista: %s', dirE);
end

%% ======================================================================
function [T, EVt, VAt] = procesar_sesion(rutaS, sid, pac, mes, cond, cfg)
S = load(rutaS);
fams = cfg.eventos.familias;
filas = {}; ev = {}; va = {};
nombres = {};
for j = 1:numel(S.X)
    X = double(S.X{j});
    N = size(X, 2);
    [~, arch] = fileparts(S.archivos{j});
    nivel = S.niveles(j);
    if strcmp(cond, 'FIS')
        E = eventos_fis(N, S.stims{j});
        % Validación del detector ERD donde la verdad es conocida
        Ee = detectar_erd(X, cfg);
        enTarea = false(1, N);
        for k = 1:height(E), enTarea(E.inicio(k):E.fin(k)) = true; end
        t0 = round(cfg.eventos.reposo_inicial_s * cfg.fs) + 1;
        frac = mean(enTarea(t0:end));
        pe = NaN; if height(Ee) > 0, pe = mean(enTarea(Ee.minimo)); end
        detect = arrayfun(@(k) any(Ee.minimo >= E.inicio(k) & Ee.minimo <= E.fin(k)), 1:height(E));
        va(end+1, :) = {sid, arch, nivel, height(Ee), pe, frac, mean(detect)}; %#ok<AGROW>
    else
        E = detectar_erd(X, cfg);
    end
    ev(end+1, :) = {sid, arch, cond, nivel, height(E)}; %#ok<AGROW>
    malas = zeros(0, 2);
    if isfield(S, 'malas'), malas = S.malas{j}; end
    for def = {'A', 'B'}
        W = ventanas_por_definicion(E, N, def{1}, cond, cfg);
        L = W.L(1);
        W = W(ventana_valida(malas, W.inicio, L), :);      % rechazo de artefactos
        if isempty(W), continue, end
        % referencia Δ: media del reposo inicial del archivo, misma longitud L
        sr = ventanas_reposo_inicial(L, cfg);
        sr = sr(ventana_valida(malas, sr, L));
        Fr = [];
        for i = 1:numel(sr)
            Fr(i, :) = caracterizadores_limites(X(:, sr(i):sr(i)+L-1), fams, cfg); %#ok<AGROW>
        end
        if isempty(Fr), ref = NaN; else, ref = mean(Fr, 1, 'omitnan'); end
        for i = 1:height(W)
            [f, nombres] = caracterizadores_limites(X(:, W.inicio(i):W.inicio(i)+L-1), fams, cfg);
            filas(end+1, :) = {sid, pac, mes, cond, arch, nivel, def{1}, L, W.tipo{i}, W.evento(i), ...
                W.inicio(i), f, f - ref}; %#ok<AGROW>
        end
    end
end
T = cell2table(filas, 'VariableNames', {'sesion','paciente','mes','condicion','archivo','nivel_archivo', ...
    'definicion','L','tipo','evento','inicio','feat','feat_delta'});
EVt = cell2table(ev, 'VariableNames', {'sesion','archivo','condicion','nivel_archivo','n_eventos'});
if isempty(va)
    VAt = table();
else
    VAt = cell2table(va, 'VariableNames', {'sesion','archivo','nivel_archivo','n_erd', ...
        'prop_eventos_en_tarea','frac_tiempo_tarea','prop_bloques_detectados'});
end
end

function U = construir_datasets(F, nombres, cfg)
% Escribe los ARFF de cada prueba y devuelve las unidades de evaluación.
U = struct('clave', {}, 'arff', {}, 'meta', {}, 'info', {});
pruebas = {'3c', '2c', 'huella', 'tarea_rep'};
for def = {'A', 'B'}
    for var = {'crudo', 'delta'}
        for cond = cfg.experimentos.condiciones
            Fc = F(strcmp(F.definicion, def{1}) & strcmp(F.condicion, cond{1}), :);
            if isempty(Fc), continue, end
            for p = pruebas
                [Fp, clases] = etiquetar(Fc, p{1});
                if isempty(Fp), continue, end
                sesiones = unique(Fp.sesion);
                bal = cell(numel(sesiones), 1);
                for s = 1:numel(sesiones)
                    bal{s} = balancear(Fp(strcmp(Fp.sesion, sesiones{s}), :));
                end
                alcances = [sesiones; {'GLOBAL'}];
                for a = 1:numel(alcances)
                    if a <= numel(sesiones)
                        D = bal{a};
                        if isempty(D) || height(D) < cfg.eventos.min_por_clase * numel(clases), continue, end
                    else
                        D = vertcat(bal{:});
                        if isempty(D), continue, end
                    end
                    clave = sprintf('%s_%s_%s_%s_%s', def{1}, var{1}, p{1}, cond{1}, alcances{a});
                    dirA = fullfile(cfg.rutas.resultados, 'arff', sprintf('reposo%s_%s_%s', ...
                        cfg.eventos.sufijo, def{1}, var{1}), p{1});
                    ruta = fullfile(dirA, sprintf('%s_%s.arff', cond{1}, alcances{a}));
                    if strcmp(var{1}, 'crudo'), X = a_matriz(D.feat); else, X = a_matriz(D.feat_delta); end
                    if ~isfile(ruta)
                        exportar_arff(ruta, clave, nombres, X, clases, D.etiqueta, ...
                            {sprintf('Etapa 8: definicion %s, variante %s, prueba %s, L = %d muestras', ...
                            def{1}, var{1}, p{1}, D.L(1))});
                        writetable(D(:, {'sesion','archivo','nivel_archivo','tipo','evento','inicio','etiqueta'}), ...
                            strrep(ruta, '.arff', '_meta.csv'));
                    end
                    info = struct('definicion', def{1}, 'variante', var{1}, 'prueba', p{1}, ...
                        'condicion', cond{1}, 'alcance', alcances{a}, 'L', D.L(1));
                    U(end+1) = struct('clave', clave, 'arff', ruta, 'meta', strrep(ruta, '.arff', '_meta.csv'), 'info', info); %#ok<AGROW>
                end
            end
        end
    end
end
end

function X = a_matriz(c)
% La columna de características puede llegar como celda de filas o como matriz.
if iscell(c), X = vertcat(c{:}); else, X = c; end
end

function [D, clases] = etiquetar(F, prueba)
D = F;
esRep = strcmp(D.tipo, 'reposo');
switch prueba
    case '3c'
        clases = {'N00','N10','N40'};
        D.etiqueta = repmat({'N00'}, height(D), 1);
        D.etiqueta(~esRep) = arrayfun(@(n) sprintf('N%02d', n), D.nivel_archivo(~esRep), 'UniformOutput', false);
    case '2c'
        clases = {'N10','N40'};
        D = D(~esRep, :);
        D.etiqueta = arrayfun(@(n) sprintf('N%02d', n), D.nivel_archivo, 'UniformOutput', false);
    case 'huella'
        clases = {'R10','R40'};
        D = D(esRep, :);
        D.etiqueta = arrayfun(@(n) sprintf('R%02d', n), D.nivel_archivo, 'UniformOutput', false);
    case 'tarea_rep'
        clases = {'R','T'};
        D.etiqueta = repmat({'R'}, height(D), 1);
        D.etiqueta(~esRep) = {'T'};
end
end

function D = balancear(D)
% Submuestreo determinista (equiespaciado) de cada clase al tamaño de la menor.
if isempty(D), return, end
[g, ~] = findgroups(D.etiqueta);
n = min(accumarray(g, 1));
if numel(unique(g)) < 2 || n == 0, D = D([], :); return, end
mantener = false(height(D), 1);
for k = 1:max(g)
    idx = find(g == k);
    [~, o] = sortrows(D(idx, {'archivo','evento','inicio'}));
    idx = idx(o);
    mantener(idx(unique(round(linspace(1, numel(idx), n))))) = true;
end
D = D(mantener, :);
end

function evaluar_unidad(u, modelos, etiq, cfg, ruta_cache)
if isfile(ruta_cache), return, end
weka_iniciar(cfg);
data = weka_cargar_arff(u.arff);
meta = readtable(u.meta, 'TextType', 'char');
conjuntos = cfg.eventos.conjuntos;
if ~(strcmp(u.info.variante, 'crudo') && strcmp(u.info.prueba, '3c'))
    conjuntos = conjuntos(1);                 % la referencia del artículo solo en 3c crudo
end
filas = {};
for c = 1:numel(conjuntos)
    fam = strsplit(conjuntos{c}, '+');
    attrs = {};
    for f = 1:numel(fam), attrs = [attrs, strcat(fam{f}, '_', cfg.canales.nombres)]; end %#ok<AGROW>
    [D, nest] = weka_seleccionar_atributos(data, attrs);
    if nest == 0, continue, end
    for i = 1:numel(modelos)
        try
            [pred, real, clases] = weka_predicciones_cv(D, modelos(i), cfg, cfg.eventos.semillas);
        catch err
            filas(end+1, :) = {u.info.definicion, u.info.variante, u.info.prueba, u.info.condicion, ...
                u.info.alcance, u.info.L, conjuntos{c}, etiq{i}, 0, NaN, NaN, NaN, NaN, NaN, NaN, ['error: ' err.message]}; %#ok<AGROW>
            continue
        end
        acc = 100 * mean(pred == real, 1);
        nc = numel(clases);
        C = zeros(nc);
        for j = 1:size(pred, 2), C = C + accumarray([real, pred(:, j)], 1, [nc nc]); end
        tpr = diag(C).' ./ sum(C, 2).';
        herencia = NaN; tpr_rep = NaN;
        if strcmp(u.info.prueba, '3c')
            tpr_rep = tpr(1);
            rep = strcmp(meta.tipo, 'reposo');
            P = pred(rep, :); niv = repmat(meta.nivel_archivo(rep), 1, size(pred, 2));
            comoTarea = P > 1;
            coincide = (P == 2 & niv == 10) | (P == 3 & niv == 40);
            herencia = 100 * nnz(coincide) / max(nnz(comoTarea), 1);
        end
        filas(end+1, :) = {u.info.definicion, u.info.variante, u.info.prueba, u.info.condicion, ...
            u.info.alcance, u.info.L, conjuntos{c}, etiq{i}, numel(real), mean(acc), std(acc), ...
            100 * mean(tpr), 100 * tpr_rep, herencia, 100 / nc, mat2str(C)}; %#ok<AGROW>
    end
end
if isempty(filas), filas = cell(0, 16); end            % dataset sin conjuntos estimables
R = cell2table(filas, 'VariableNames', {'definicion','variante','prueba','condicion','alcance','L', ...
    'conjunto','modelo','n','acc','acc_de_semillas','acc_balanceada','tpr_reposo','herencia_reposo', ...
    'azar','confusion'}); %#ok<NASGU>
tmp = [ruta_cache '.tmp.mat']; save(tmp, 'R'); movefile(tmp, ruta_cache, 'f');
end

function T = resumen_final(R)
R = R(~startsWith(R.confusion, 'error'), :);
intra = R(~strcmp(R.alcance, 'GLOBAL'), :);
[g, d, v, p, c, cj, m] = findgroups(intra.definicion, intra.variante, intra.prueba, intra.condicion, ...
    intra.conjunto, intra.modelo);
T = table(d, v, p, c, cj, m, splitapply(@numel, intra.acc, g), splitapply(@mean, intra.acc, g), ...
    splitapply(@std, intra.acc, g), splitapply(@mean, intra.acc_balanceada, g), ...
    splitapply(@(x) mean(x, 'omitnan'), intra.tpr_reposo, g), ...
    splitapply(@(x) mean(x, 'omitnan'), intra.herencia_reposo, g), splitapply(@mean, intra.azar, g), ...
    'VariableNames', {'definicion','variante','prueba','condicion','conjunto','modelo','n_sesiones', ...
    'acc_intra','acc_intra_de','acc_bal_intra','tpr_reposo','herencia_reposo','azar'});
G = R(strcmp(R.alcance, 'GLOBAL'), :);
T.acc_global = nan(height(T), 1);
for i = 1:height(T)
    k = find(strcmp(G.definicion, T.definicion{i}) & strcmp(G.variante, T.variante{i}) & ...
        strcmp(G.prueba, T.prueba{i}) & strcmp(G.condicion, T.condicion{i}) & ...
        strcmp(G.conjunto, T.conjunto{i}) & strcmp(G.modelo, T.modelo{i}), 1);
    if ~isempty(k), T.acc_global(i) = G.acc(k); end
end
T = sortrows(T, {'prueba','condicion','definicion','variante','acc_intra'}, {'ascend','ascend','ascend','ascend','descend'});
end

function graficar_eventos(T, dirEF)
pr = {'huella', 'tarea_rep', '2c', '3c'};
f = figure('Visible', 'off', 'Position', [100 100 1300 420]);
for i = 1:numel(pr)
    subplot(1, numel(pr), i);
    s = T(strcmp(T.prueba, pr{i}) & strcmp(T.conjunto, 'Ro64a2+Re128a2') & contains(T.modelo, 'ajustado'), :);
    if isempty(s), continue, end
    [g, d, v, c] = findgroups(s.definicion, s.variante, s.condicion);
    etiquetas = strcat(c, '-', d, '-', v);
    bar(categorical(etiquetas), splitapply(@max, s.acc_intra, g)); hold on
    yline(mean(s.azar), '--', 'azar'); grid on
    ylabel('Mejor precisión intra (%)'); title(pr{i}, 'Interpreter', 'none');
end
exportgraphics(f, fullfile(dirEF, 'pruebas_reposo.png'), 'Resolution', 150);
close(f);
end
