function acc = loso_por_sesion(T, nombres, m, cfg, evaluar)
%LOSO_POR_SESION  LOSO 10 vs 40 con una fila por sesión de T (NaN si no es evaluable).
%   acc = loso_por_sesion(T, nombres, m, cfg, evaluar). Igual que
%   loso_niveles (z-score por sesión, entrena con el resto de sesiones),
%   pero devuelve numel(unique(T.sesion, 'stable')) valores, alineados con
%   esas sesiones, para poder cruzarlos con otros criterios por sesión.
%   evaluar (opcional): sesiones que se evalúan (las demás quedan NaN).
ses = unique(T.sesion, 'stable');
if nargin < 5 || isempty(evaluar), evaluar = ses; end
acc = nan(numel(ses), 1);
X = double(matriz_caracteristicas(T));
for s = 1:numel(ses)
    i = strcmp(T.sesion, ses{s});
    X(i, :) = (X(i, :) - mean(X(i, :), 1, 'omitnan')) ./ std(X(i, :), 0, 1, 'omitnan');
end
y = arrayfun(@(n) sprintf('N%02d', n), T.nivel, 'UniformOutput', false);
weka_iniciar(cfg);
for s = 1:numel(ses)
    te = strcmp(T.sesion, ses{s});
    if numel(unique(y(te))) < 2 || ~ismember(ses{s}, evaluar), continue, end
    Dtr = weka_instancias('tr', nombres, X(~te, :), {'N10','N40'}, y(~te));
    Dte = weka_instancias('te', nombres, X(te, :), {'N10','N40'}, y(te));
    cls = weka_clasificador_param(m, cfg); cls.buildClassifier(Dtr);
    ok = 0;
    for i = 0:Dte.numInstances()-1
        ok = ok + (cls.classifyInstance(Dte.instance(i)) == Dte.instance(i).classValue());
    end
    acc(s) = 100 * ok / Dte.numInstances();
end
end
