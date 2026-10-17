%% COMPARAR_VARIANTES  Puentes x rechazo de artefactos: ¿qué combinación generaliza mejor?
%
%   Factor 1 - electrodos puenteados (r > 0.99 tras CAR):
%     P0 sin corrección, todas las sesiones (corrida CAR actual)
%     P1 interpolar puentes          + excluir sesiones con > 6 puenteados
%     P2 eliminar puentes (NaN)      + excluir sesiones con > 6 puenteados
%     P3 sin corrección              + excluir sesiones con > 6 puenteados
%   Factor 2 - rechazo de artefactos (sobre la señal limpia CAR + ICA):
%     R0 sin rechazo | R1 6 MAD, 1 canal | R2 8 MAD, 1 canal |
%     R3 6 MAD, >= 2 canales | R4 8 MAD, >= 2 canales
%   Métricas (las que NO puede inflar la huella de cada archivo van primero):
%     1) dejar-una-sesión-fuera (LOSO), 10 vs 40 %, KMI y FIS: precisión
%        media, sesiones > 50 %, Wilcoxon vs 50 %
%     2) Px010-M1-KMI (orden invertido) en LOSO
%     3) excedente intra-sesión: (10 vs 40 %) - (control de mitades)
%     4) % de ventanas conservadas
%   Ventanas contiguas de 512 muestras en todo el registro; caracterizadores
%   calibrados con CAR (Re64a2 + Re128a4); SVM ajustado; z-score por sesión.
%   Paso final: prueba de permutación (etiquetas 10/40 intercambiadas al azar
%   por sesión completa) para la mejor combinación y la actual (P0-R1).
%   -> resultados_var/*.csv

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_car');
cfg = config_niveles();
dirV = fullfile(raiz, 'resultados_var');
if ~isfolder(dirV), mkdir(dirV); end
umbral_sesion = 6;
fams = struct('nombre', {'Re64a2', 'Re128a4'}, 'modo', {'ecuacion13', 'ecuacion13'}, ...
    'm_max', {64, 128}, 'm_min', {2, 4});
rech = struct('id', {'R0','R1','R2','R3','R4'}, 'umbral_z', {Inf, 6, 8, 6, 8}, 'min_canales', {1, 1, 1, 2, 2});
puentes = struct('id', {'P0','P1','P2','P3'}, 'metodo', {'ninguno','interpolar','eliminar','ninguno'}, ...
    'excluir', {false, true, true, true});
nw = usar_paralelo(cfg);
registrar(cfg, '==== Comparación de variantes: %d puentes x %d rechazos ====', numel(puentes), numel(rech));

%% 1) Preprocesamiento por método de puentes (P1, P2); P0/P3 reutilizan resultados_car
dirsSen = struct('ninguno', fullfile(cfg.rutas.resultados, 'senales'));
for m = {'interpolar', 'eliminar'}
    c = cfg;
    c.rutas.resultados = fullfile(dirV, ['puentes_' m{1}]);
    c.puentes.metodo = m{1};
    if ~isfolder(c.rutas.resultados), mkdir(c.rutas.resultados); end
    etapa1_preprocesar(c);
    dirsSen.(m{1}) = fullfile(c.rutas.resultados, 'senales');
end

%% 2) Sesiones admitidas por la regla de calidad
Pp = readtable(fullfile(cfg.rutas.resultados, 'puentes_por_archivo.csv'), 'TextType', 'char');
Pp.sesion = strcat(Pp.paciente, '_M', string(Pp.mes), '_', Pp.condicion);
g = groupsummary(Pp, 'sesion', 'max', 'n_puenteados');
excluidas = g.sesion(g.max_n_puenteados > umbral_sesion);
registrar(cfg, '  sesiones excluidas por > %d canales puenteados: %s', umbral_sesion, strjoin(excluidas, ', '));
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
ses = ses(logical(ses.completa), :);

%% 3) Ventanas y características por método de puentes (con validez por rechazo)
W = struct();
for m = fieldnames(dirsSen).'
    arch = fullfile(dirV, sprintf('ventanas_%s.mat', m{1}));
    if isfile(arch), S = load(arch); W.(m{1}) = S.T; continue, end
    partes = cell(height(ses), 1);
    sid = ses.sesion; cnd = ses.condicion; dS = dirsSen.(m{1});
    parfor (s = 1:height(ses), nw)
        partes{s} = ventanas_tiempo_sesion(fullfile(dS, [sid{s} '.mat']), sid{s}, cnd{s}, fams, 512, cfg, rech);
    end
    T = vertcat(partes{:});
    save(arch, 'T', '-v7.3');
    W.(m{1}) = T;
    registrar(cfg, '  ventanas (%s): %d', m{1}, height(T));
end

%% 4) Evaluación de las 20 combinaciones
H = readtable(fullfile(raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), 'Delimiter', ',', 'TextType', 'char');
H.Properties.VariableNames = {'modelo','especificacion'};
modelo = modelo_desde_id(H.especificacion{strcmp(H.modelo, 'SVM (ajustado)')});
nombres = {};
for k = 1:numel(fams), nombres = [nombres, strcat(fams(k).nombre, '_', cfg.canales.nombres)]; end %#ok<AGROW>
combos = {};
for p = 1:numel(puentes), for r = 1:numel(rech), combos(end+1, :) = {p, r}; end, end %#ok<AGROW>
res = cell(size(combos, 1), 1);
parfor (k = 1:size(combos, 1), nw)
    p = combos{k, 1}; r = combos{k, 2};
    T = W.(puentes(p).metodo);
    v = T.valida; if iscell(v), v = vertcat(v{:}); end
    T = T(v(:, r), :);
    if puentes(p).excluir, T = T(~ismember(T.sesion, excluidas), :); end
    res{k} = evaluar_combinacion(T, nombres, modelo, cfg, puentes(p).id, rech(r).id, ...
        height(T) / height(W.(puentes(p).metodo)));
end
R = vertcat(res{:});
R = sortrows(R, {'loso_kmi','excedente_kmi'}, {'descend','descend'});
writetable(R, fullfile(dirV, 'comparacion_variantes.csv'));
disp(R);

%% 5) Prueba de permutación: mejor combinación vs actual (P0-R1)
mejor = R(1, :);
cand = unique([{[mejor.puentes{1} '-' mejor.rechazo{1}]}, {'P0-R1'}], 'stable');
perm = {};
for c = cand
    pr = strsplit(c{1}, '-');
    p = find(strcmp({puentes.id}, pr{1})); r = find(strcmp({rech.id}, pr{2}));
    T = W.(puentes(p).metodo); v = T.valida; if iscell(v), v = vertcat(v{:}); end; T = T(v(:, r), :);
    if puentes(p).excluir, T = T(~ismember(T.sesion, excluidas), :); end
    T = T(strcmp(T.condicion, 'KMI'), :);
    obs = loso_niveles(T, nombres, modelo, cfg, []);
    np = 100; nula = zeros(np, 1);
    sesK = unique(T.sesion);
    parfor (i = 1:np, nw)
        s = RandStream('mt19937ar', 'Seed', 500 + i);
        voltear = sesK(rand(s, numel(sesK), 1) < 0.5);       % sesiones con etiquetas 10/40 intercambiadas
        nula(i) = mean(loso_niveles(T, nombres, modelo, cfg, voltear));
    end
    pval = (1 + nnz(nula >= mean(obs))) / (1 + np);
    perm(end+1, :) = {c{1}, mean(obs), mean(nula), std(nula), percentil95(nula), pval}; %#ok<AGROW>
    registrar(cfg, '  permutación %s: LOSO KMI = %.1f %% | nula %.1f +- %.1f (p95 %.1f) | p = %.3f', ...
        c{1}, mean(obs), mean(nula), std(nula), percentil95(nula), pval);
end
P = cell2table(perm, 'VariableNames', {'combinacion','loso_kmi','nula_media','nula_de','nula_p95','p'});
writetable(P, fullfile(dirV, 'permutacion_loso.csv'));
registrar(cfg, '==== Comparación de variantes terminada ====');

%% ======================================================================
function fila = evaluar_combinacion(T, nombres, modelo, cfg, pid, rid, frac)
fila = [table({pid}, {rid}, 100 * frac, 'VariableNames', {'puentes','rechazo','pct_ventanas'}), ...
    evaluar_generalizacion(T, nombres, modelo, cfg)];
end
