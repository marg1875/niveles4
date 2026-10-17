function [Mm, Ms] = weka_cv_repetida(data, m, cfg, semillas)
%WEKA_CV_REPETIDA  CV k-fold repetida con varias semillas.
%   [Mm, Ms] = weka_cv_repetida(data, m, cfg, semillas) evalúa la
%   especificación m (o un nombre de clasificador) con una CV k-fold por
%   cada semilla (reutiliza weka_validar) y devuelve la media (Mm) y la
%   desviación estándar entre repeticiones (Ms) de acc, kappa, mae, wf, y de
%   las métricas por clase. Mm.confusion es la suma de las matrices.
%   Con semillas = 1 coincide exactamente con weka_validar.
campos = {'acc','kappa','mae','rmse','wtpr','wfpr','wprec','wrec','wf','wroc'};
porclase = {'tpr','fpr','prec','rec','f','roc'};
R = cell(numel(semillas), 1);
for i = 1:numel(semillas)
    c = cfg; c.weka.semilla = semillas(i); c.weka.cv_modo = 'estratificada';
    if isstruct(m), cls = weka_clasificador_param(m, c); else, cls = m; end
    R{i} = weka_validar(data, cls, c);
end
R = [R{:}];
for f = [campos, porclase]
    X = vertcat(R.(f{1}));
    Mm.(f{1}) = mean(X, 1);
    Ms.(f{1}) = std(X, 0, 1);
end
Mm.clases = R(1).clases;
Mm.n = R(1).n;
Mm.confusion = sum(cat(3, R.confusion), 3);
Mm.acc_repeticiones = [R.acc];
end
