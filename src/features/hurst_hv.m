function [H, nvalidos] = hurst_hv(V, cfg)
%HURST_HV  Hurst con semivariograma (HV, Ec. 16 + R/S).
%   H = hurst_hv(V, cfg), V (C x N). Igual que el script del grupo
%   (Hurst_Variograma_Serie_Tiempo_PM10.m, l. 396 y 479):
%       [vario, alt] = variograma(v, tam);  [tau, RS] = rangors(vario, 9);
%   es decir, se calcula gamma(h), h = 1..floor(N/2), y se aplica el R/S
%   con particiones (p = cfg.hurst.hv_particiones) a la secuencia gamma.
%   OJO: Niveles3 calculaba HV como pendiente/2 del log-variograma, lo que
%   NO corresponde al artículo.
N = size(V, 2);
G = semivariograma(V, floor(N/2));
[H, nvalidos] = hurst_hrs(G, cfg.hurst.hv_particiones, cfg);
end
