function acc = loso_niveles(T, nombres, modelo, cfg, voltear, evaluar)
%LOSO_NIVELES  Exactitud 10 vs 40 % dejando una sesión fuera (Weka).
%   acc = loso_niveles(T, nombres, modelo, cfg, voltear, evaluar)
%   T: tabla de ventanas (ventanas_tiempo_sesion) con sesion, nivel, feat.
%   Cada sesión se normaliza con su propio z-score (sin usar etiquetas), se
%   entrena con el resto de sesiones y se evalúa en ella. Devuelve la
%   exactitud (%) por sesión evaluable (las que tienen ambos niveles).
%   'voltear': sesiones cuyas etiquetas 10/40 se intercambian (prueba de
%   permutación por sesión completa); [] o {} = etiquetas reales.
%   evaluar (opcional): sesiones de prueba (por defecto todas); el
%   entrenamiento usa siempre el resto de sesiones de T.
X = matriz_caracteristicas(T);
ses = unique(T.sesion, 'stable');
if nargin < 6 || isempty(evaluar), evaluar = ses; end
for s = 1:numel(ses)
    i = strcmp(T.sesion, ses{s});
    X(i, :) = (X(i, :) - mean(X(i, :), 1, 'omitnan')) ./ std(X(i, :), 0, 1, 'omitnan');
end
y = T.nivel;
if isempty(voltear), voltear = {}; end
v = ismember(T.sesion, voltear);
y(v) = 50 - y(v);                                  % 10 <-> 40
y = arrayfun(@(n) sprintf('N%02d', n), y, 'UniformOutput', false);
weka_iniciar(cfg);
acc = nan(numel(ses), 1);
for s = 1:numel(ses)
    te = strcmp(T.sesion, ses{s});
    if numel(unique(y(te))) < 2 || ~ismember(ses{s}, evaluar), continue, end
    Dtr = weka_instancias('tr', nombres, X(~te, :), {'N10','N40'}, y(~te));   % en memoria
    Dte = weka_instancias('te', nombres, X(te, :), {'N10','N40'}, y(te));
    cls = weka_clasificador_param(modelo, cfg); cls.buildClassifier(Dtr);
    ok = 0;
    for i = 0:Dte.numInstances()-1
        ok = ok + (cls.classifyInstance(Dte.instance(i)) == Dte.instance(i).classValue());
    end
    acc(s) = 100 * ok / Dte.numInstances();
end
acc = acc(isfinite(acc));
end
