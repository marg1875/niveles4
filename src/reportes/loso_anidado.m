function [acc, elegidos] = loso_anidado(T, nombres, conjuntos, modelo, cfg, nw)
%LOSO_ANIDADO  LOSO 10 vs 40 % con la ELECCIÓN del conjunto dentro de cada fold.
%   [acc, elegidos] = loso_anidado(T, nombres, conjuntos, modelo, cfg, nw)
%   T: ventanas (ventanas_tiempo_sesion) con todas las columnas 'nombres'
%   en T.feat. conjuntos: celda de 'FamA+FamB' (atributos '<fam>_<canal>').
%   Para cada sesión de prueba s: con las demás sesiones se calcula la LOSO
%   interna de cada conjunto, se elige el mejor y se evalúa en s
%   entrenando con todas las demás. Así la exactitud no está inflada por
%   haber elegido el conjunto con los mismos datos de prueba.
%   acc (%) y conjunto elegido por sesión evaluable. nw: trabajadores del
%   parfor sobre las sesiones de prueba (0 = en serie).
if nargin < 6, nw = 0; end
ses = unique(T.sesion, 'stable');
X = matriz_caracteristicas(T);
pref = cellfun(@(n) n(1:find(n == '_', 1, 'last') - 1), nombres, 'UniformOutput', false);
cols = cellfun(@(c) find(ismember(pref, strsplit(c, '+'))), conjuntos, 'UniformOutput', false);
acc = nan(numel(ses), 1);
elegidos = cell(numel(ses), 1);
parfor (s = 1:numel(ses), nw)
    te = strcmp(T.sesion, ses{s});
    if numel(unique(T.nivel(te))) < 2, continue, end
    Ttr = T(~te, :);
    interna = zeros(numel(conjuntos), 1);
    for c = 1:numel(conjuntos)
        Tc = Ttr; Tc.feat = X(~te, cols{c});
        interna(c) = mean(loso_niveles(Tc, nombres(cols{c}), modelo, cfg, []));
    end
    [~, j] = max(interna);
    elegidos{s} = conjuntos{j};
    Tc = T; Tc.feat = X(:, cols{j});
    acc(s) = loso_niveles(Tc, nombres(cols{j}), modelo, cfg, [], ses(s));
end
ok = isfinite(acc);
acc = acc(ok); elegidos = elegidos(ok);
end
