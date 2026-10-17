function [H, nvalidos] = hurst_hrs(V, p, cfg)
%HURST_HRS  Hurst con particiones (HRS, Ecs. 7-9 y 14-15 del artículo).
%   H = hurst_hrs(V, p, cfg), V (C x N). Usa los niveles diádicos de
%   rangors(v, p): m = N, N/2, ..., hasta el primer m < N/p.
N = size(V, 2);
ms = niveles_rangors(N, p, cfg.hurst.bloque_min_abs);
[H, nvalidos] = pendiente_hurst(ms, rs_diadico(V, ms, cfg.hurst.promedio_rs), cfg.hurst.min_escalas);
end
