function [H, nvalidos] = hurst_ho(V, cfg)
%HURST_HO  Hurst "original" (HO, Ecs. 10-15 del artículo).
%   H = hurst_ho(V, cfg), V (C x N). R/S clásico con TODAS las particiones
%   diádicas m = N, N/2, N/4, ... con bloque >= cfg.hurst.ho_bloque_min
%   (interpretación documentada de "uses all possible partitions"; por
%   debajo de ~8 muestras el estadístico R/S es muy sesgado).
N = size(V, 2);
ms = floor(N ./ 2.^(0:floor(log2(N))));
ms = ms(ms >= max(cfg.hurst.ho_bloque_min, cfg.hurst.bloque_min_abs));
[H, nvalidos] = pendiente_hurst(ms, rs_diadico(V, ms, cfg.hurst.promedio_rs), cfg.hurst.min_escalas);
end
