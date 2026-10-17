function [idx, curt] = limpiar_artefactos_ica(Y, cfg)
%LIMPIAR_ARTEFACTOS_ICA  Selecciona componentes ICA artefactuales por curtosis.
%   [idx, curt] = limpiar_artefactos_ica(Y, cfg) devuelve los índices de las
%   fuentes Y (k x n) con curtosis atípicamente alta (como máximo
%   cfg.ica.max_removidos, de mayor a menor). Los parpadeos, los movimientos
%   y las ráfagas musculares producen fuentes fuertemente super-gaussianas.
%
%   Criterio robusto: z = (log(kurt+3) - mediana) / (1.4826 * MAD) >
%   cfg.ica.curtosis_nsigma. Se usa escala logarítmica y estadísticos
%   robustos porque en estos registros la curtosis llega a valores de
%   10^3-10^4, y con media/desviación estándar una sola componente extrema
%   enmascara a las demás.
%
%   NOTA: las versiones previas (ica_.m, Niveles3/ica.py) eliminaban la
%   componente de "mayor varianza" tras blanquear; en ese espacio todas las
%   varianzas valen ~1, por lo que esa elección era arbitraria.
curt = mean(Y.^4, 2) ./ mean(Y.^2, 2).^2 - 3;
lk = log(max(curt + 3, eps));
mad = median(abs(lk - median(lk)));
z = (lk - median(lk)) ./ max(1.4826 * mad, eps);
cand = find(z > cfg.ica.curtosis_nsigma);
[~, orden] = sort(curt(cand), 'descend');
idx = cand(orden(1:min(numel(cand), cfg.ica.max_removidos)));
idx = idx(:).';
end
