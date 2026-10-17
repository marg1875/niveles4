function [p, chi2, gl] = prueba_friedman(A)
%PRUEBA_FRIEDMAN  Prueba de Friedman (bloques = filas, tratamientos = columnas).
%   [p, chi2, gl] = prueba_friedman(A), A (n sesiones x k configuraciones).
%   Se eliminan las filas con NaN. Implementación propia (sin Statistics
%   Toolbox): rangos por fila con empates promediados, estadístico con
%   corrección por empates y p a partir de chi-cuadrado con k-1 gl
%   (gammainc de MATLAB base).
A = A(all(isfinite(A), 2), :);
[n, k] = size(A);
gl = k - 1;
if n < 2 || k < 2, p = NaN; chi2 = NaN; return, end
Rg = zeros(n, k);
for i = 1:n, Rg(i,:) = rangos_promedio(A(i,:)); end
Rj = sum(Rg, 1);
num = 12/(n*k*(k+1)) * sum(Rj.^2) - 3*n*(k+1);
% corrección por empates
T = 0;
for i = 1:n
    [~, ~, g] = unique(A(i,:));
    t = accumarray(g(:), 1);
    T = T + sum(t.^3 - t);
end
den = 1 - T / (n*(k^3 - k));
chi2 = num / max(den, eps);
p = gammainc(chi2/2, gl/2, 'upper');
end

function r = rangos_promedio(x)
[~, o] = sort(x);
r = zeros(size(x));
r(o) = 1:numel(x);
[u, ~, g] = unique(x);
for j = 1:numel(u)
    sel = g == j;
    r(sel) = mean(r(sel));
end
end
