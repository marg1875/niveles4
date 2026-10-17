%% PRUEBA_ORDEN_INVERTIDO  ¿El modelo aprende nivel o tiempo? (Px010-M1, KMI)
%
%   En casi todas las sesiones el orden fue FIS10 -> KMI10 -> FIS40 -> KMI40:
%   el 40 % siempre va después. En Px010-M1 la imaginación al 40 % se grabó
%   ANTES que la del 10 %. Entrenando con todas las demás sesiones KMI
%   (dejar una sesión fuera) y probando en Px010-M1:
%     - si el modelo aprendió TIEMPO ("lo posterior = 40 %"), invertirá las
%       etiquetas (precisión << 50 %);
%     - si aprendió NIVEL, acertará (precisión > 50 %).
%   Como referencia se hace lo mismo con cada sesión KMI de orden normal.
%   Ventanas: las de analisis_deriva_temporal.m (512 muestras, todo el
%   registro). Características: conjunto calibrado. Variantes: crudo y
%   estandarizado por sesión (z-score de cada característica con TODAS las
%   ventanas de la sesión, sin usar etiquetas).
%   -> resultados/deriva/prueba_orden_invertido.csv

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
cfg = config_niveles();
weka_iniciar(cfg);
dirD = fullfile(cfg.rutas.resultados, 'deriva');
S = load(fullfile(dirD, 'ventanas_tiempo.mat'));
W = S.W; nombres = S.nombres;
if iscell(W.feat), X = vertcat(W.feat{:}); else, X = W.feat; end
conjunto = cfg.eventos.conjuntos{1};
cols = startsWith(nombres, strcat(strsplit(conjunto, '+'), '_'));
X = X(:, cols); nom = nombres(cols);
k = strcmp(W.condicion, 'KMI');
W = W(k, :); X = X(k, :);

% Orden temporal de cada sesión: ¿el 40 % se grabó antes que el 10 %?
ses = unique(W.sesion, 'stable');
orden = cell(numel(ses), 1);
for s = 1:numel(ses)
    i = strcmp(W.sesion, ses{s});
    t10 = min(W.t_min(i & W.nivel == 10)); t40 = min(W.t_min(i & W.nivel == 40));
    if t40 < t10, orden{s} = '40 antes que 10 (invertido)'; else, orden{s} = '10 antes que 40 (normal)'; end
end

H = readtable(fullfile(cfg.rutas.resultados, 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), 'Delimiter', ',', 'TextType', 'char');
H.Properties.VariableNames = {'modelo','especificacion'};
modelos = [modelo_desde_id(H.especificacion{strcmp(H.modelo, 'SVM (ajustado)')}), ...
           modelo_desde_id(H.especificacion{strcmp(H.modelo, 'MLP (ajustado)')}), ...
           especificar_modelo('SVM', struct())];
etiq = {'SVM ajustado', 'MLP ajustado', 'SVM artículo'};

filas = {};
for variante = {'crudo', 'z_por_sesion'}
    Xv = X;
    if strcmp(variante{1}, 'z_por_sesion')
        for s = 1:numel(ses)
            i = strcmp(W.sesion, ses{s});
            Xv(i, :) = (X(i, :) - mean(X(i, :), 1, 'omitnan')) ./ std(X(i, :), 0, 1, 'omitnan');
        end
    end
    y = arrayfun(@(n) sprintf('N%02d', n), W.nivel, 'UniformOutput', false);
    for s = 1:numel(ses)
        te = strcmp(W.sesion, ses{s}); tr = ~te;
        rtr = [tempname '.arff']; rte = [tempname '.arff'];
        exportar_arff(rtr, 'entrena', nom, Xv(tr, :), {'N10','N40'}, y(tr), {});
        exportar_arff(rte, 'prueba', nom, Xv(te, :), {'N10','N40'}, y(te), {});
        Dtr = weka_cargar_arff(rtr); Dte = weka_cargar_arff(rte);
        for m = 1:numel(modelos)
            cls = weka_clasificador_param(modelos(m), cfg);
            cls.buildClassifier(Dtr);
            p = zeros(Dte.numInstances(), 1);
            for i = 0:Dte.numInstances()-1, p(i+1) = cls.classifyInstance(Dte.instance(i)); end
            real = strcmp(y(te), 'N40');
            pred40 = p == 1;
            acc = 100 * mean(pred40 == real);
            f10 = 100 * mean(pred40(~real)); f40 = 100 * mean(pred40(real));
            filas(end+1, :) = {ses{s}, orden{s}, variante{1}, etiq{m}, acc, f10, f40}; %#ok<AGROW>
        end
        delete(rtr); delete(rte);
    end
end
T = cell2table(filas, 'VariableNames', {'sesion_prueba','orden','variante','modelo','acc', ...
    'pct_pred40_en_archivo10','pct_pred40_en_archivo40'});
writetable(T, fullfile(dirD, 'prueba_orden_invertido.csv'));

fprintf('\n=== Sesión de orden INVERTIDO (Px010-M1) ===\n');
disp(T(contains(T.orden, 'invertido'), :));
fprintf('\n=== Sesiones de orden normal (media de dejar-una-fuera) ===\n');
N = T(~contains(T.orden, 'invertido'), :);
disp(groupsummary(N, {'variante','modelo'}, {'mean','std','min','max'}, 'acc'));
