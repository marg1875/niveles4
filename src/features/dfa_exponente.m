function A = dfa_exponente(V, cfg, ns, orden)
%DFA_EXPONENTE  Exponente alfa del análisis de fluctuaciones sin tendencia.
%   A = dfa_exponente(V, cfg), V (C x N). Port de compute_dfa() (Niveles3;
%   Peng et al., 1995): escalas n = unique(floor(logspace(log10(n_min),
%   log10(0.33 N), 15))), n_min = max(16, 0.064 fs), con >= 3 ventanas;
%   DFA1 (destendencia lineal, dfa_curva); alfa = pendiente de log10 F vs
%   log10 n (>= 5 escalas, si no NaN).
%   A = dfa_exponente(V, cfg, ns, orden): escalas y orden explícitos
%   (calibración, etapa 10).
%   Ruido blanco -> 0.5; ruido 1/f -> 1; movimiento browniano -> 1.5.
p = cfg.fractales;
[C, N] = size(V);
A = nan(C, 1);
if nargin < 3
    nmin = max(floor(p.dfa_escala_min_s * cfg.fs), 16);
    nmax = floor(N * p.dfa_escala_max_frac);
    if nmax <= nmin, return, end
    ns = unique(floor(logspace(log10(nmin), log10(nmax), p.dfa_n_escalas)));
    ns = ns(ns > 0 & floor(N ./ ns) >= p.dfa_min_ventanas);
    orden = 1;
end
ok = std(double(V), 1, 2) >= 1e-8;
F = dfa_curva(V, ns, orden);
A(ok) = pendiente_hurst(ns, F(ok, :), p.dfa_min_escalas);
end
