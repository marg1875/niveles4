function H = variograma_pendiente(V, cfg)
%VARIOGRAMA_PENDIENTE  Hurst por la pendiente del semivariograma (VAR de Niveles3).
%   H = variograma_pendiente(V, cfg), V (C x N). Port de compute_variogram():
%     gamma(h) = 1/(2(N-h)) sum (x(t+h) - x(t))^2 (semivariograma.m),
%     h = 1 .. h_max - 1, h_max = min(max(floor(0.2 N), 10), floor(N/5));
%     H = pendiente(log10 gamma vs log10 h) / 2 (>= 5 puntos, si no NaN).
%   Es distinto de HV (hurst_hv.m, script del grupo), que aplica R/S con
%   particiones a la secuencia gamma(h).
p = cfg.fractales;
N = size(V, 2);
hmax = min(max(floor(N * p.var_lag_max_frac), p.var_min_lags), floor(N / 5));
h = 1:hmax-1;
if isempty(h), H = nan(size(V, 1), 1); return, end
G = semivariograma(V, hmax - 1);
G(G <= 1e-10) = NaN;
H = pendiente_hurst(h, G, p.var_min_puntos) / 2;
end
