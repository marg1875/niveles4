function D = dimensiones_directas(V)
%DIMENSIONES_DIRECTAS  Dimensiones fractales y descriptores sin escalas, por canal.
%   D = dimensiones_directas(V), V (C x N) -> D (C x 7):
%     1 KATZ   Katz (1988): log10(n) / (log10(n) + log10(d/Lc)), Lc = suma |dx|,
%              d = max |x - x(1)|, n = N - 1.
%     2 PETd   Petrosian (1995) con cambios de signo de la derivada.
%     3 PETm   Petrosian con cruces de la mediana.
%     4 SEV    Sevcik (1998): gráfica en el cuadrado unitario, 1 + ln(Lc)/ln(2(N-1)).
%     5 HJMOB  movilidad de Hjorth sqrt(var(dx)/var(x)).
%     6 HJCOM  complejidad de Hjorth movilidad(dx)/movilidad(x).
%     7 LL     line length normalizada: media |dx| / std(x).
X = double(V);
N = size(X, 2);
dx = diff(X, 1, 2); ddx = diff(dx, 1, 2);
n = N - 1;
Lc = sum(abs(dx), 2);
d = max(abs(X - X(:, 1)), [], 2);
D = nan(size(X, 1), 7);
D(:, 1) = log10(n) ./ (log10(n) + log10(max(d, eps) ./ max(Lc, eps)));
nd = sum(diff(sign(dx), 1, 2) ~= 0, 2);
D(:, 2) = log10(N) ./ (log10(N) + log10(N ./ (N + 0.4 * nd)));
s = sign(X - median(X, 2)); nm = sum(diff(s, 1, 2) ~= 0, 2);
D(:, 3) = log10(N) ./ (log10(N) + log10(N ./ (N + 0.4 * nm)));
Xn = (X - min(X, [], 2)) ./ max(max(X, [], 2) - min(X, [], 2), eps);
Ls = sum(sqrt(diff(Xn, 1, 2).^2 + (1 / n)^2), 2);
D(:, 4) = 1 + log(Ls) / log(2 * n);
vx = var(X, 0, 2); vd = var(dx, 0, 2); vdd = var(ddx, 0, 2);
D(:, 5) = sqrt(vd ./ vx);
D(:, 6) = sqrt(vdd ./ vd) ./ D(:, 5);
D(:, 7) = mean(abs(dx), 2) ./ std(X, 0, 2);
D(~isfinite(D)) = NaN;
end
