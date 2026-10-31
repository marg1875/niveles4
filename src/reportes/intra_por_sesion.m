function [acc, K] = intra_por_sesion(T, nombres, m, cfg, semillas, k)
%INTRA_POR_SESION  Exactitud 10 vs 40 dentro de cada sesión (CV estratificada).
%   [acc, K] = intra_por_sesion(T, nombres, m, cfg, semillas, k). Para cada
%   sesión de T, validación cruzada estratificada de k particiones sobre sus
%   ventanas, repetida con cada semilla; acc(s) es la media sobre semillas.
%   Mismo z-score por sesión y mismo clasificador m que loso_por_sesion,
%   para que ambas cifras sean comparables. NaN si la sesión no tiene ambas
%   clases. K(:, :, s): matriz de confusión de la sesión s (fila real 10/40,
%   columna predicha 10/40), sumada sobre las semillas.
ses = unique(T.sesion, 'stable');
acc = nan(numel(ses), 1);
K = zeros(2, 2, numel(ses));
X = double(matriz_caracteristicas(T));
y = arrayfun(@(n) sprintf('N%02d', n), T.nivel, 'UniformOutput', false);
weka_iniciar(cfg);
for s = 1:numel(ses)
    i = find(strcmp(T.sesion, ses{s}));
    if numel(unique(y(i))) < 2, continue, end
    Xs = X(i, :);
    Xs = (Xs - mean(Xs, 1, 'omitnan')) ./ std(Xs, 0, 1, 'omitnan');
    ys = y(i);
    a = nan(numel(semillas), 1);
    for r = 1:numel(semillas)
        rng(semillas(r), 'twister');
        cv = cvpartition(ys, 'KFold', k);
        ok = 0;
        for f = 1:k
            tr = training(cv, f); te = test(cv, f);
            Dtr = weka_instancias('tr', nombres, Xs(tr, :), {'N10','N40'}, ys(tr));
            Dte = weka_instancias('te', nombres, Xs(te, :), {'N10','N40'}, ys(te));
            cls = weka_clasificador_param(m, cfg); cls.buildClassifier(Dtr);
            for q = 0:Dte.numInstances()-1
                yr = Dte.instance(q).classValue(); yp = cls.classifyInstance(Dte.instance(q));
                ok = ok + (yp == yr);
                K(yr + 1, yp + 1, s) = K(yr + 1, yp + 1, s) + 1;
            end
        end
        a(r) = 100 * ok / numel(ys);
    end
    acc(s) = mean(a);
end
end
