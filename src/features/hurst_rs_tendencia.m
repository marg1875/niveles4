function H = hurst_rs_tendencia(V, cfg)
%HURST_RS_TENDENCIA  Hurst R/S con eliminación de tendencia (RS de Niveles3).
%   H = hurst_rs_tendencia(V, cfg), V (C x N). Port de compute_rs()
%   (Niveles3/src/features/fractal.py):
%     1) se resta la media y la recta de mínimos cuadrados de cada canal;
%     2) escalas n = unique(round(logspace(log10(n_min), log10(0.40 N), 20))),
%        n_min = max(32, 0.128 fs); se exigen >= 5 bloques por escala;
%     3) <R/S>_n = media sobre bloques de R/S, con S = desviación estándar
%        poblacional (ddof = 0) y bloques constantes omitidos;
%     4) H = pendiente de log10<R/S> vs log10 n (>= 5 escalas, si no NaN).
%   Diferencias con HO/HRS (rs_diadico): escalas logarítmicas y no
%   diádicas, destendencia previa y S poblacional.
p = cfg.fractales;
[C, N] = size(V);
H = nan(C, 1);
if N < 40, return, end
X = double(V);
t = (0:N-1) - (N-1)/2;
X = X - mean(X, 2);
X = X - (X * t.' / sum(t.^2)) .* t;                % destendencia lineal
ok = std(X, 1, 2) >= 1e-10;
nmin = max(floor(p.rs_escala_min_s * cfg.fs), 32);
nmax = floor(N * p.rs_escala_max_frac);
if nmax <= nmin, return, end
ns = unique(round(logspace(log10(nmin), log10(nmax), p.rs_n_escalas)));
ns = ns(ns >= 4 & ns < N & floor(N ./ ns) >= p.rs_min_bloques);
RS = rs_bloques(X, ns);
H(ok) = pendiente_hurst(ns, RS(ok, :), p.rs_min_escalas);
end

function RS = rs_bloques(X, ns)
% <R/S>_n por canal (bloques contiguos, S poblacional, bloques planos fuera).
[C, N] = size(X);
RS = nan(C, numel(ns));
Xt = X.';
for i = 1:numel(ns)
    n = ns(i); nb = floor(N / n);
    B = reshape(Xt(1:nb*n, :), n, nb, C);
    S = std(B, 1, 1);
    Y = cumsum(B - mean(B, 1), 1);
    R = max(Y, [], 1) - min(Y, [], 1);
    q = R ./ S;
    q(S < 1e-10) = NaN;
    RS(:, i) = squeeze(mean(q, 2, 'omitnan'));
end
end
