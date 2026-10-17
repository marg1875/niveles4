function D = higuchi_dimension(V, cfg, kmin, kmax)
%HIGUCHI_DIMENSION  Dimensión fractal de Higuchi (1988) por canal.
%   D = higuchi_dimension(V, cfg), V (C x N). Port de compute_higuchi()
%   (Niveles3): k_max = min(floor(0.25 N), 30) con piso 5, k = 1..k_max;
%   D = -pendiente de log10 L(k) vs log10 k (higuchi_curva; >= 5 puntos
%   válidos, si no NaN).
%   D = higuchi_dimension(V, cfg, kmin, kmax): rango de k explícito
%   (calibración, etapa 10).
%   Ruido blanco -> D ~ 2; movimiento browniano -> D ~ 1.5; curva suave -> D ~ 1.
p = cfg.fractales;
N = size(V, 2);
if nargin < 4
    kmin = 1;
    kmax = max(min(floor(N * p.hig_kmax_frac), p.hig_kmax_tope), p.hig_kmax_piso);
end
Lk = higuchi_curva(V, kmax);
D = -pendiente_hurst(kmin:kmax, Lk(:, kmin:kmax), p.hig_min_puntos);
end
