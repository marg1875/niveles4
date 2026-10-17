function M = rejillas_hiperparametros(cfg, clf)
%REJILLAS_HIPERPARAMETROS  Modelos a evaluar en la búsqueda de la etapa 6.
%   M = rejillas_hiperparametros(cfg) devuelve un arreglo de especificaciones
%   (especificar_modelo) para kNN, SVM, MLP y RF, cada uno cruzado con las
%   opciones de selección de características cfg.ajuste.selecciones.
%   M = rejillas_hiperparametros(cfg, clf) solo las de ese clasificador.
%   Las rejillas se editan en cfg.ajuste.rejilla (config_niveles.m).
if nargin < 2, clf = {'kNN','SVM','MLP','RF'}; end
clf = cellstr(clf);
R = cfg.ajuste.rejilla;
P = {};
if ismember('kNN', clf)
    for K = R.knn_K, for d = R.knn_dist, for w = R.knn_peso
        P{end+1} = {'kNN', struct('K', K, 'dist', d{1}, 'peso', w{1})}; %#ok<AGROW>
    end, end, end
end
if ismember('SVM', clf)
    for f = R.svm_filtro
        for E = R.svm_E, for C = R.svm_C_poly
            P{end+1} = {'SVM', struct('kernel', 'poly', 'E', E, 'C', C, 'filtro', f)}; %#ok<AGROW>
        end, end
        for g = R.svm_gamma, for C = R.svm_C_rbf
            P{end+1} = {'SVM', struct('kernel', 'rbf', 'gamma', g, 'C', C, 'filtro', f)}; %#ok<AGROW>
        end, end
    end
end
if ismember('MLP', clf)
    for H = R.mlp_H, for L = R.mlp_L, for N = R.mlp_N, for V = R.mlp_V
        P{end+1} = {'MLP', struct('H', H{1}, 'L', L, 'N', N, 'V', V)}; %#ok<AGROW>
    end, end, end, end
end
if ismember('RF', clf)
    for a = R.rf_arboles, for K = R.rf_K, for pr = R.rf_prof
        P{end+1} = {'RF', struct('arboles', a, 'K', K, 'prof', pr)}; %#ok<AGROW>
    end, end, end
end
M = [];
for i = 1:numel(P)
    for s = cfg.ajuste.selecciones
        M = [M, especificar_modelo(P{i}{1}, P{i}{2}, s{1})]; %#ok<AGROW>
    end
end
end
