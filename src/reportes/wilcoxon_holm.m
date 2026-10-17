function T = wilcoxon_holm(A, nombres)
%WILCOXON_HOLM  Comparaciones por pares (Wilcoxon de rangos con signo) + Holm.
%   T = wilcoxon_holm(A, nombres), A (n sesiones x k), nombres (1 x k).
%   Devuelve una tabla con cada par, la mediana de la diferencia, p crudo y
%   p ajustado por Holm-Bonferroni. Usa signrank (Statistics Toolbox) si
%   está disponible; si no, la aproximación normal con corrección por empates.
k = size(A, 2);
pares = nchoosek(1:k, 2);
np = size(pares, 1);
a = cell(np,1); b = a; med = nan(np,1); p = nan(np,1);
usar_signrank = exist('signrank', 'file') == 2;
for i = 1:np
    x = A(:, pares(i,1)); y = A(:, pares(i,2));
    ok = isfinite(x) & isfinite(y);
    d = x(ok) - y(ok);
    a{i} = nombres{pares(i,1)}; b{i} = nombres{pares(i,2)};
    med(i) = median(d);
    if nnz(d ~= 0) < 1, p(i) = 1; continue, end
    if usar_signrank
        p(i) = signrank(x(ok), y(ok));
    else
        p(i) = signrank_normal(d);
    end
end
% Holm
[ps, o] = sort(p);
padj = zeros(size(p));
acum = 0;
for j = 1:np
    acum = max(acum, min(1, (np - j + 1) * ps(j)));
    padj(o(j)) = acum;
end
T = table(a, b, med, p, padj, 'VariableNames', {'config_a','config_b','mediana_dif','p','p_holm'});
end

function p = signrank_normal(d)
d = d(d ~= 0);
n = numel(d);
[~, o] = sort(abs(d));
r = zeros(n,1); r(o) = 1:n;
[u, ~, g] = unique(abs(d));
for j = 1:numel(u), r(g==j) = mean(r(g==j)); end
W = sum(r(d > 0));
mu = n*(n+1)/4;
t = accumarray(g, 1);
s = sqrt(n*(n+1)*(2*n+1)/24 - sum(t.^3 - t)/48);
z = (W - mu - 0.5*sign(W - mu)) / s;
p = erfc(abs(z)/sqrt(2));
end
