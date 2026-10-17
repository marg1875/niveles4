%% COMPARAR_ASR  Reparación de picos con ASR frente a la señal sin reparar.
%   Análisis principal SIN rechazo de ventanas y con todas las sesiones:
%     A0      CAR + ICA (resultados_car), sin ASR
%     ASR20   ASR k = 20 antes del CAR (resultados_var/asr_k20)
%     ASR30   ASR k = 30 (resultados_var/asr_k30)
%   y como sensibilidad A0 con el rechazo de 6 MAD / 1 canal (A0-R1).
%   Para cada variante: generalización entre sesiones 10 vs 40 % (LOSO, SVM
%   ajustado, Re64a2 + Re128a4, ventanas de 512), Wilcoxon vs 50 %, Px010-M1
%   invertido, excedente sobre el control de mitades, % de ventanas que el
%   rechazo de 6 MAD seguiría marcando (eficacia de la reparación) y % del
%   registro reparado por ASR.
%   Prueba de permutación final (1000 permutaciones de etiquetas 10/40 por
%   sesión completa) en KMI y FIS para A0 y la mejor variante ASR.
%   -> resultados_var/comparacion_asr.csv, resultados_var/permutacion_final.csv

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_car');
cfg = config_niveles();
dirV = fullfile(raiz, 'resultados_var');
if ~isfolder(dirV), mkdir(dirV); end
n_perm = 1000;
lote = 50;                                         % permutaciones por lote (avance)
ks = [20 30];
fams = struct('nombre', {'Re64a2', 'Re128a4'}, 'modo', {'ecuacion13', 'ecuacion13'}, ...
    'm_max', {64, 128}, 'm_min', {2, 4});
rech = struct('id', {'R0','R1'}, 'umbral_z', {Inf, 6}, 'min_canales', {1, 1});
nw = usar_paralelo(cfg);
registrar(cfg, '==== Comparación ASR: sin ASR vs k = %s ====', mat2str(ks));

%% 1) Preprocesamiento con ASR (sin corrección de puentes, como A0)
dirsSen = struct('A0', fullfile(cfg.rutas.resultados, 'senales'));
for k = ks
    c = cfg;
    c.rutas.resultados = fullfile(dirV, sprintf('asr_k%d', k));
    c.asr.activar = true;
    c.asr.k = k;
    if ~isfolder(c.rutas.resultados), mkdir(c.rutas.resultados); end
    etapa1_preprocesar(c);
    dirsSen.(sprintf('ASR%d', k)) = fullfile(c.rutas.resultados, 'senales');
end
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
ses = ses(logical(ses.completa), :);

%% 2) Ventanas (A0 reutiliza las de comparar_variantes; columna 2 de 'valida' = R1)
W = struct();
pctASR = struct('A0', 0);
for m = fieldnames(dirsSen).'
    if strcmp(m{1}, 'A0')
        arch = fullfile(dirV, 'ventanas_ninguno.mat');
    else
        arch = fullfile(dirV, sprintf('ventanas_%s.mat', lower(m{1})));
    end
    if isfile(arch)
        S = load(arch); T = S.T;
    else
        partes = cell(height(ses), 1);
        sid = ses.sesion; cnd = ses.condicion; dS = dirsSen.(m{1});
        parfor (s = 1:height(ses), nw)
            partes{s} = ventanas_tiempo_sesion(fullfile(dS, [sid{s} '.mat']), sid{s}, cnd{s}, fams, 512, cfg, rech);
        end
        T = vertcat(partes{:});
        save(arch, 'T', '-v7.3');
    end
    v = T.valida; if iscell(v), v = vertcat(v{:}); end
    T.valida = v(:, 1:2);                          % [R0 R1]
    W.(m{1}) = T;
    if ~strcmp(m{1}, 'A0')
        fr = cellfun(@(s) mean(getfield(load(fullfile(dirsSen.(m{1}), [s '.mat']), 'info'), 'info', 'frac_asr')), ses.sesion); %#ok<GFLD>
        pctASR.(m{1}) = 100 * mean(fr);
    end
    registrar(cfg, '  ventanas %s: %d (R1 marcaría %.1f %%)', m{1}, height(T), 100 * mean(~T.valida(:, 2)));
end

%% 3) Evaluación de las variantes
H = readtable(fullfile(raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), 'Delimiter', ',', 'TextType', 'char');
H.Properties.VariableNames = {'modelo','especificacion'};
modelo = modelo_desde_id(H.especificacion{strcmp(H.modelo, 'SVM (ajustado)')});
nombres = {};
for k = 1:numel(fams), nombres = [nombres, strcat(fams(k).nombre, '_', cfg.canales.nombres)]; end %#ok<AGROW>
variantes = [fieldnames(dirsSen); {'A0-R1'}];
res = cell(numel(variantes), 1);
parfor (i = 1:numel(variantes), nw)
    base = strtok(variantes{i}, '-');
    T = W.(base);
    marc = 100 * mean(~T.valida(:, 2));
    if contains(variantes{i}, '-R1'), T = T(T.valida(:, 2), :); end
    res{i} = [table(variantes(i), 100 * height(T) / height(W.(base)), pctASR.(base), marc, ...
        'VariableNames', {'variante','pct_ventanas','pct_reparado_asr','pct_marcaria_R1'}), ...
        evaluar_generalizacion(T, nombres, modelo, cfg)];
end
R = vertcat(res{:});
writetable(R, fullfile(dirV, 'comparacion_asr.csv'));
disp(R);
for i = 1:height(R)
    registrar(cfg, '  %-6s LOSO KMI %.1f %% (%d/%d > 50, p=%.3f) | FIS %.1f %% | Px010 %.0f %% | ventanas %.0f %%', ...
        R.variante{i}, R.loso_kmi(i), R.ses_sobre50_kmi(i), R.n_ses_kmi(i), R.p_wilcoxon_kmi(i), ...
        R.loso_fis(i), R.px010_invertido(i), R.pct_ventanas(i));
end

%% 4) Permutación final: A0 y la mejor variante ASR, KMI y FIS
Ra = R(startsWith(R.variante, 'ASR'), :);
[~, im] = max(Ra.loso_kmi);
cand = {'A0', Ra.variante{im}};
perm = {};
for c = cand
    for cond = {'KMI', 'FIS'}
        T = W.(c{1});
        T = T(strcmp(T.condicion, cond{1}), :);
        obs = mean(loso_niveles(T, nombres, modelo, cfg, []));
        sesC = unique(T.sesion);
        nula = zeros(n_perm, 1);
        fav = fullfile(dirV, sprintf('permutacion_%s_%s.avance.txt', c{1}, cond{1}));
        for b = 1:lote:n_perm
            idx = b:min(n_perm, b + lote - 1);
            parte = zeros(numel(idx), 1);
            parfor (j = 1:numel(idx), nw)
                s = RandStream('mt19937ar', 'Seed', 500 + idx(j));
                voltear = sesC(rand(s, numel(sesC), 1) < 0.5);   % sesiones con 10/40 intercambiados
                parte(j) = mean(loso_niveles(T, nombres, modelo, cfg, voltear));
            end
            nula(idx) = parte;
            fid = fopen(fav, 'w'); fprintf(fid, '%d %d', idx(end), n_perm); fclose(fid);
        end
        pval = (1 + nnz(nula >= obs)) / (1 + n_perm);
        perm(end+1, :) = {c{1}, cond{1}, obs, mean(nula), std(nula), percentil95(nula), pval}; %#ok<AGROW>
        registrar(cfg, '  permutación final %s %s: LOSO = %.1f %% | nula %.1f +- %.1f (p95 %.1f) | p = %.4f', ...
            c{1}, cond{1}, obs, mean(nula), std(nula), percentil95(nula), pval);
    end
end
P = cell2table(perm, 'VariableNames', {'variante','condicion','loso','nula_media','nula_de','nula_p95','p'});
writetable(P, fullfile(dirV, 'permutacion_final.csv'));
registrar(cfg, '==== Comparación ASR terminada ====');
