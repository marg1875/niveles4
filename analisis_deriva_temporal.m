%% ANALISIS_DERIVA_TEMPORAL  ¿La separación 10 % vs 40 % es nivel o tiempo (fatiga)?
%
%   Las grabaciones de cada sesión se hicieron seguidas, casi siempre en el
%   orden FIS10 -> KMI10 -> FIS40 -> KMI40 (~4 min entre archivos, según la
%   marca "Created on" del encabezado de cada .mat). El 40 % siempre va
%   ~8 min después que el 10 %: nivel y tiempo/fatiga están confundidos.
%
%   Dentro de un archivo el nivel es constante y el tiempo avanza, así que
%   la pendiente de cada caracterizador vs tiempo DENTRO de los archivos
%   (beta, por minuto; modelo y = a_archivo + beta t) mide la deriva
%   temporal sin confundido de nivel. Con ella:
%     (a) diferencia predicha por el tiempo = beta * dt (dt = separación
%         entre los centros de los dos archivos) frente a la diferencia
%         observada D = media(40 %) - media(10 %);
%     (b) control sin nivel: clasificar la 1.a vs la 2.a mitad del MISMO
%         archivo (solo tiempo, ~1.3 min) y compararlo con 10 vs 40 %
%         (~8 min). Si el control es tan alto como 10 vs 40 %, la
%         no estacionariedad temporal basta para explicar esa precisión
%         con CV por ventanas.
%   Resultado de (a): la pendiente intra-archivo NO extrapola a la escala de
%   la sesión (no es una fatiga lineal), por eso no se resta.
%   Ventanas: no traslapadas de 512 muestras en todo el registro (desde
%   0.5 s). Caracterizadores: los calibrados (Ro64a2, Re128a2) y HO/HRS.
%   Nivel por archivo: si los dos archivos de la sesión tienen stims
%   distintos (10 y 40) se usa stims (corrige Px008-M3-KMI, cuyos nombres
%   están intercambiados); si no, el nombre del archivo.
%   -> resultados/deriva/*.csv y figuras

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
cfg = config_niveles();
warning('off', 'MATLAB:nearlySingularMatrix');   % sesiones con archivos reguardados (hora atípica)
dirD = fullfile(cfg.rutas.resultados, 'deriva');
if ~isfolder(dirD), mkdir(dirD); end
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
ses = ses(logical(ses.completa), :);
fams = cfg.eventos.familias;
L = 512;
nw = usar_paralelo(cfg);

%% 1) Características por ventana con su tiempo absoluto en la sesión
partes = cell(height(ses), 1);
sid = ses.sesion; cond = ses.condicion;
parfor (s = 1:height(ses), nw)
    partes{s} = ventanas_tiempo(fullfile(cfg.rutas.resultados, 'senales', [sid{s} '.mat']), ...
        sid{s}, cond{s}, fams, L, cfg);
end
W = vertcat(partes{:});
nombres = W.Properties.UserData;
if isempty(nombres)
    for k = 1:numel(fams), nombres = [nombres, strcat(fams(k).nombre, '_', cfg.canales.nombres)]; end %#ok<AGROW>
end
save(fullfile(dirD, 'ventanas_tiempo.mat'), 'W', 'nombres', '-v7.3');

%% 2) Deriva dentro de archivos vs diferencia entre archivos
if iscell(W.feat), X = vertcat(W.feat{:}); else, X = W.feat; end
filas = {};
for s = 1:height(ses)
    idx = strcmp(W.sesion, ses.sesion{s});
    Ws = W(idx, :); Xs = X(idx, :);
    i10 = Ws.nivel == 10; i40 = Ws.nivel == 40;
    if ~any(i10) || ~any(i40), continue, end
    dt = mean(Ws.t_min(i40)) - mean(Ws.t_min(i10));
    for j = 1:size(Xs, 2)
        y = Xs(:, j);
        ok = isfinite(y);
        if nnz(ok & i10) < 5 || nnz(ok & i40) < 5, continue, end
        % y = a_10 * [10] + a_40 * [40] + beta * t   (pendiente común DENTRO de archivos)
        A = [double(i10(ok)), double(i40(ok)), Ws.t_min(ok)];
        [b, se] = mco(A, y(ok));
        D = mean(y(ok & i40)) - mean(y(ok & i10));
        sD = sqrt(var(y(ok & i40))/nnz(ok & i40) + var(y(ok & i10))/nnz(ok & i10));
        pred = b(3) * dt;
        z = (D - pred) / sqrt(sD^2 + (se(3) * dt)^2);
        filas(end+1, :) = {ses.sesion{s}, ses.condicion{s}, nombres{j}, dt, D, b(3), pred, D - pred, z}; %#ok<AGROW>
    end
end
T = cell2table(filas, 'VariableNames', {'sesion','condicion','caracteristica','dt_min', ...
    'dif_observada','pendiente_por_min','dif_predicha_tiempo','dif_residual_nivel','z_residual'});
writetable(T, fullfile(dirD, 'deriva_por_caracteristica.csv'));

% Resumen: ¿qué fracción de la diferencia observada explica el tiempo?
T.fam = strtok(T.caracteristica, '_');
R = groupsummary(T, {'condicion', 'fam'}, {'median'}, {'dif_observada', 'dif_predicha_tiempo', 'z_residual'});
T.misma_direccion = sign(T.dif_observada) == sign(T.dif_predicha_tiempo);
T.frac_explicada = min(abs(T.dif_predicha_tiempo) ./ abs(T.dif_observada), 5) .* T.misma_direccion;
T.residual_signif = abs(T.z_residual) > 1.96;
R2 = groupsummary(T, {'condicion', 'fam'}, {'mean'}, {'misma_direccion', 'frac_explicada', 'residual_signif'});
Rres = join(R, R2(:, {'condicion','fam','mean_misma_direccion','mean_frac_explicada','mean_residual_signif'}));
writetable(Rres, fullfile(dirD, 'deriva_resumen.csv'));
disp(Rres);

%% 3) Clasificación 10 vs 40 %: cruda vs sin tendencia temporal
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), 'Delimiter', ',', 'TextType', 'char');
H.Properties.VariableNames = {'modelo','especificacion'};
modelos = [modelo_desde_id(H.especificacion{strcmp(H.modelo, 'SVM (ajustado)')}), ...
           modelo_desde_id(H.especificacion{strcmp(H.modelo, 'MLP (ajustado)')}), ...
           especificar_modelo('SVM', struct())];
etiq = {'SVM ajustado', 'MLP ajustado', 'SVM artículo'};
conjunto = cfg.eventos.conjuntos{1};
selCols = startsWith(nombres, strcat(strsplit(conjunto, '+'), '_'));
dirTmp = fullfile(dirD, 'arff'); if ~isfolder(dirTmp), mkdir(dirTmp); end
filas = cell(height(ses), 1);
parfor (s = 1:height(ses), nw)
    filas{s} = clasificar_sesion(W(strcmp(W.sesion, ses.sesion{s}), :), selCols, nombres, modelos, etiq, ...
        cfg, fullfile(dirTmp, ses.sesion{s}));
end
C = vertcat(filas{:});
writetable(C, fullfile(dirD, 'clasificacion_cruda_vs_destendenciada.csv'));
Cs = groupsummary(C, {'condicion', 'modelo', 'variante'}, {'mean', 'std'}, 'acc');
writetable(Cs, fullfile(dirD, 'clasificacion_resumen.csv'));
disp(Cs);
graficar_deriva(W, X, nombres, T, C, dirD);
fprintf('Resultados en %s\n', dirD);

%% ======================================================================
function T = ventanas_tiempo(rutaS, sid, cond, fams, L, cfg)
S = load(rutaS);
fs = cfg.fs;
% Hora de cada archivo (encabezado del .mat = momento en que se guardó,
% al final de la grabación) y nivel corregido
hora = NaT(1, numel(S.X)); nivel = S.niveles;
stimsMax = cellfun(@max, S.stims);
if numel(S.X) == 2 && all(sort(stimsMax) == [10 40])
    nivel = stimsMax;                         % stims distintos y coherentes: se confía en stims
end
for j = 1:numel(S.X)
    hora(j) = leer_hora_mat(S.archivos{j});
end
t0 = min(hora - seconds(cellfun(@(x) size(x, 2), S.X) / fs));   % inicio de la primera grabación
filas = {};
for j = 1:numel(S.X)
    X = double(S.X{j});
    N = size(X, 2);
    inicio_archivo = hora(j) - seconds(N / fs);
    off_min = minutes(inicio_archivo - t0);
    [~, arch] = fileparts(S.archivos{j});
    starts = round(0.5 * fs) + 1 : L : N - L + 1;
    if isfield(S, 'malas'), starts = starts(ventana_valida(S.malas{j}, starts, L)); end
    for s = starts
        f = caracterizadores_limites(X(:, s:s+L-1), fams, cfg);
        t = off_min + (s + L/2) / fs / 60;
        filas(end+1, :) = {sid, cond, arch, nivel(j), (s + L/2) / fs / 60, t, f}; %#ok<AGROW>
    end
end
T = cell2table(filas, 'VariableNames', {'sesion','condicion','archivo','nivel','t_archivo_min','t_min','feat'});
end

function h = leer_hora_mat(ruta)
fid = fopen(ruta, 'r'); b = fread(fid, 116, '*char')'; fclose(fid);
tok = regexp(b, 'Created on: (\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})', 'tokens', 'once');
h = datetime(regexprep(tok{1}, ' +', ' '), 'InputFormat', 'eee MMM d HH:mm:ss yyyy', 'Locale', 'en_US');
end

function [b, se] = mco(A, y)
% Mínimos cuadrados con errores estándar.
b = A \ y;
r = y - A * b;
gl = max(numel(y) - size(A, 2), 1);
s2 = (r' * r) / gl;
se = sqrt(diag(s2 * inv(A' * A)));
end

function R = clasificar_sesion(Ws, selCols, nombres, modelos, etiq, cfg, base)
weka_iniciar(cfg);
if iscell(Ws.feat), X = vertcat(Ws.feat{:}); else, X = Ws.feat; end
X = X(:, selCols); nom = nombres(selCols);
y = arrayfun(@(n) sprintf('N%02d', n), Ws.nivel, 'UniformOutput', false);
% NOTA: una versión anterior restaba la tendencia temporal intra-archivo
% extrapolada a toda la sesión; se descartó porque esa pendiente no
% extrapola (predice diferencias ~10x mayores y de signo opuesto al
% observado) e introducía separaciones artificiales.
%
% Variantes:
%   'nivel_10_vs_40'  : ventanas del archivo de 10 % vs del de 40 % (~8 min de separación)
%   'control_mitades' : 1.a vs 2.a mitad del MISMO archivo (~1.3 min; sin nivel)
%                       -> cuánto separa la sola no estacionariedad temporal
filas = {};
archivos = unique(Ws.archivo, 'stable');
casos = {'nivel_10_vs_40', {}};
for a = 1:numel(archivos), casos(end+1, :) = {'control_mitades', archivos(a)}; end %#ok<AGROW>
for c = 1:size(casos, 1)
    if strcmp(casos{c, 1}, 'nivel_10_vs_40')
        M = X; yy = y; clases = {'N10', 'N40'};
    else
        sel = strcmp(Ws.archivo, casos{c, 2}{1});
        M = X(sel, :);
        t = Ws.t_archivo_min(sel);
        yy = repmat({'M1'}, nnz(sel), 1); yy(t > median(t)) = {'M2'};
        clases = {'M1', 'M2'};
    end
    % Con el rechazo de artefactos algunos archivos quedan con muy pocas
    % ventanas: se exige al menos cfg.weka.folds por clase.
    if min(cellfun(@(c) nnz(strcmp(yy, c)), clases)) < cfg.weka.folds, continue, end
    ruta = sprintf('%s_%s_%d.arff', base, casos{c, 1}, c);
    exportar_arff(ruta, 'deriva', nom, M, clases, yy, {});
    D = weka_cargar_arff(ruta);
    for i = 1:numel(modelos)
        [pred, real] = weka_predicciones_cv(D, modelos(i), cfg, 1:5);
        filas(end+1, :) = {Ws.sesion{1}, Ws.condicion{1}, etiq{i}, casos{c, 1}, ...
            100 * mean(pred(:) == repmat(real, 5, 1))}; %#ok<AGROW>
    end
end
if isempty(filas)          % sesión sin casos evaluables (pocas ventanas tras el rechazo)
    R = table(cell(0,1), cell(0,1), cell(0,1), cell(0,1), zeros(0,1), ...
        'VariableNames', {'sesion','condicion','modelo','variante','acc'});
else
    R = cell2table(filas, 'VariableNames', {'sesion','condicion','modelo','variante','acc'});
end
end

function graficar_deriva(W, X, nombres, T, C, dirD)
f = figure('Visible', 'off', 'Position', [100 100 1300 420]);
% (1) Ejemplo: media de Ro64a2 (todos los canales) vs tiempo en la sesión
subplot(1, 3, 1); hold on
cols = startsWith(nombres, 'Ro64a2_');
for s = unique(W.sesion(strcmp(W.condicion, 'KMI'))).'
    i = strcmp(W.sesion, s{1});
    scatter(W.t_min(i), mean(X(i, cols), 2, 'omitnan'), 6, W.nivel(i), 'filled');
end
colormap([0 0.45 0.74; 0.85 0.33 0.1]); grid on
xlabel('Minutos desde el inicio de la sesión'); ylabel('Ro64a2 medio');
title('KMI: caracterizador vs tiempo (azul 10 %, rojo 40 %)');
% (2) diferencia observada vs predicha por el tiempo
subplot(1, 3, 2);
t = T(startsWith(T.caracteristica, {'Ro64a2_', 'Re128a2_'}), :);
scatter(t.dif_predicha_tiempo, t.dif_observada, 8, strcmp(t.condicion, 'KMI'), 'filled'); hold on
lim = max(abs([t.dif_predicha_tiempo; t.dif_observada])); plot([-lim lim], [-lim lim], 'k--'); grid on
xlabel('Diferencia 40-10 predicha solo por el tiempo'); ylabel('Diferencia 40-10 observada');
title('Puntos sobre la diagonal = el tiempo lo explica');
% (3) precisión cruda vs sin tendencia
subplot(1, 3, 3);
[g, c, m, v] = findgroups(C.condicion, C.modelo, C.variante);
acc = splitapply(@mean, C.acc, g);
bar(categorical(strcat(c, " ", m, " ", v)), acc); yline(50, '--', 'azar'); grid on
ylabel('Precisión 10 vs 40 % (intra)'); title('Clasificación cruda vs sin tendencia temporal');
exportgraphics(f, fullfile(dirD, 'deriva_temporal.png'), 'Resolution', 150);
close(f);
end
