function fila = evaluar_generalizacion(T, nombres, modelo, cfg)
%EVALUAR_GENERALIZACION  Métricas 10 vs 40 % entre y dentro de sesiones.
%   fila = evaluar_generalizacion(T, nombres, modelo, cfg), T tabla de
%   ventanas (ventanas_tiempo_sesion). Para KMI y FIS:
%     loso_*          exactitud media dejando una sesión fuera (loso_niveles)
%     n_ses_*, ses_sobre50_*, p_wilcoxon_* (vs 50 %)
%     px010_invertido (KMI; sesión con orden de registro 40 -> 10)
%     excedente_*     CV intra-sesión de nivel menos el control de mitades
%                     del mismo archivo (deriva temporal).
fila = table();
for cond = {'KMI', 'FIS'}
    c = lower(cond{1});
    Tc = T(strcmp(T.condicion, cond{1}), :);
    acc = loso_niveles(Tc, nombres, modelo, cfg, []);
    ses = unique(Tc.sesion, 'stable');
    fila.(['loso_' c]) = mean(acc);
    fila.(['n_ses_' c]) = numel(acc);
    fila.(['ses_sobre50_' c]) = nnz(acc > 50);
    W = wilcoxon_holm([acc(:), 50 * ones(numel(acc), 1)], {'a', 'b'});
    fila.(['p_wilcoxon_' c]) = W.p;
    if strcmp(cond{1}, 'KMI')
        i = find(strcmp(ses, 'Px010_M1_KMI'));
        if isempty(i), fila.px010_invertido = NaN; else, fila.px010_invertido = acc(i); end
    end
    ex = nan(numel(ses), 1);
    for s = 1:numel(ses)
        Ts = Tc(strcmp(Tc.sesion, ses{s}), :);
        a = intra_cv(Ts, nombres, modelo, cfg, 'nivel');
        b = mean(cellfun(@(f) intra_cv(Ts(strcmp(Ts.archivo, f), :), nombres, modelo, cfg, 'mitad'), unique(Ts.archivo)), 'omitnan');
        ex(s) = a - b;
    end
    fila.(['excedente_' c]) = mean(ex, 'omitnan');
end
end

function a = intra_cv(T, nombres, modelo, cfg, tipo)
% CV 5-fold dentro de una sesión (nivel) o de un archivo (mitad 1 vs 2).
if strcmp(tipo, 'nivel')
    y = arrayfun(@(n) sprintf('N%02d', n), T.nivel, 'UniformOutput', false); cl = {'N10','N40'};
else
    y = repmat({'M1'}, height(T), 1); y(T.t_archivo_min > median(T.t_archivo_min)) = {'M2'}; cl = {'M1','M2'};
end
a = NaN;
if min(cellfun(@(c) nnz(strcmp(y, c)), cl)) < cfg.weka.folds, return, end
weka_iniciar(cfg);
r = [tempname '.arff'];
exportar_arff(r, 'cv', nombres, matriz_caracteristicas(T), cl, y, {});
D = weka_cargar_arff(r);
[pred, real] = weka_predicciones_cv(D, modelo, cfg, 1);
a = 100 * mean(pred == real);
delete(r);
end
