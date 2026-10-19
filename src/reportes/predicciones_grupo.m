function pred = predicciones_grupo(T, nombres, m, cfg, grupo, voltear, entrenables, featfun)
%PREDICCIONES_GRUPO  Predicción 10/40 de cada ventana dejando fuera un grupo a la vez.
%   pred = predicciones_grupo(T, nombres, m, cfg, grupo, voltear, entrenables, featfun)
%   T: ventanas con sesion, nivel y feat (o featfun). grupo: celda con el
%   grupo de cada fila (sesión -> LOSO, paciente -> LOPO, día -> LODO).
%   Para cada grupo g se entrena con las filas 'entrenables' de los demás
%   grupos y se predicen TODAS las filas de g. Cada sesión se normaliza con
%   su propio z-score (sin etiquetas), como en loso_niveles.
%   voltear: sesiones cuyas etiquetas 10/40 se intercambian (permutación).
%   entrenables (opcional, lógico): filas que pueden entrenar (por
%   defecto todas). Así se entrena con KMI y se predice también la FIS.
%   featfun (opcional): @(entrena, y) -> matriz de características calculada
%   solo con las filas de entrenamiento y sus etiquetas y (ya permutadas),
%   p. ej. filtros CSP por fold.
%   pred: nivel predicho (10/40) por fila, NaN si no se pudo entrenar.
n = height(T);
if nargin < 6 || isempty(voltear), voltear = {}; end
if nargin < 7 || isempty(entrenables), entrenables = true(n, 1); end
if nargin < 8, featfun = []; end
y = T.nivel;
v = ismember(T.sesion, voltear);
y(v) = 50 - y(v);
ylab = arrayfun(@(k) sprintf('N%02d', k), y, 'UniformOutput', false);
ses = unique(T.sesion, 'stable');
pred = nan(n, 1);
weka_iniciar(cfg);
gs = unique(grupo, 'stable');
X0 = [];
if isempty(featfun), X0 = normalizar(double(matriz_caracteristicas(T)), T.sesion, ses); end
for i = 1:numel(gs)
    te = strcmp(grupo, gs{i});
    tr = ~te & entrenables;
    if numel(unique(y(tr))) < 2, continue, end
    if isempty(featfun), X = X0; else, X = normalizar(featfun(tr, y), T.sesion, ses); end
    Dtr = weka_instancias('tr', nombres, X(tr, :), {'N10','N40'}, ylab(tr));
    Dte = weka_instancias('te', nombres, X(te, :), {'N10','N40'}, ylab(te));
    cls = weka_clasificador_param(m, cfg); cls.buildClassifier(Dtr);
    idx = find(te);
    for j = 0:Dte.numInstances()-1
        pred(idx(j + 1)) = 10 + 30 * cls.classifyInstance(Dte.instance(j));
    end
end
end

function X = normalizar(X, sesion, ses)
for s = 1:numel(ses)
    k = strcmp(sesion, ses{s});
    X(k, :) = (X(k, :) - mean(X(k, :), 1, 'omitnan')) ./ std(X(k, :), 0, 1, 'omitnan');
end
X(~isfinite(X)) = NaN;
end
