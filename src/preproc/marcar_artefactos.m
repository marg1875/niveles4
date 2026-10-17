function malas = marcar_artefactos(X, cfg)
%MARCAR_ARTEFACTOS  Tramos con amplitud atípica en la señal ya limpia.
%   malas = marcar_artefactos(X, cfg), X (canales x muestras). Una muestra
%   es mala si en al menos cfg.artefactos.min_canales canales
%   |x - mediana| > cfg.artefactos.umbral_z * 1.4826 * MAD (desviación
%   robusta). Se extiende +-cfg.artefactos.margen_s. Los canales en NaN
%   (eliminados por puenteo) se ignoran.
%   Devuelve una matriz k x 2 con los tramos [inicio fin] (muestras, base 1).
if ~cfg.artefactos.activar
    malas = zeros(0, 2); return
end
minc = 1;
if isfield(cfg.artefactos, 'min_canales'), minc = cfg.artefactos.min_canales; end
X = X(all(isfinite(X), 2), :);
med = median(X, 2);
mad = median(abs(X - med), 2);
z = abs(X - med) ./ max(1.4826 * mad, eps);
mala = sum(z > cfg.artefactos.umbral_z, 1) >= minc;
m = round(cfg.artefactos.margen_s * cfg.fs);
if m > 0 && any(mala)
    mala = movmax(double(mala), 2*m + 1) > 0;
end
d = diff([false, mala, false]);
malas = [find(d == 1).', (find(d == -1) - 1).'];
end
