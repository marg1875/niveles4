function [X, nombres] = caracteristicas_desde_curvas(M, K, lay, cfg)
%CARACTERISTICAS_DESDE_CURVAS  Candidatos de Higuchi/DFA a partir de las curvas.
%   [X, nombres] = caracteristicas_desde_curvas(M, K, lay, cfg)
%   M (ventanas x lay.ancho): filas de curvas_ventana. K, lay: de
%   candidatos_hig_dfa. X (ventanas x numel(K)*C), nombres '<cand>_<canal>'.
%   Pendiente de log10(curva) vs log10(x) en el rango de cada candidato,
%   por mínimos cuadrados vectorizados (puntos NaN fuera; NaN si quedan
%   menos de cfg.calib_hd.min_puntos_validos).
C = lay.C; W = size(M, 1);
ch = cfg.canales.nombres;
X = nan(W, numel(K) * C);
nombres = cell(1, numel(K) * C);
curvas = struct();
for b = lay.bloques
    curvas.(b.nombre) = log10(reshape(M(:, b.o + (1:b.K*C)), [], C, b.K));
end
for k = 1:numel(K)
    cols = (k-1)*C + (1:C);
    nombres(cols) = strcat(K(k).nombre, '_', ch);
    if strcmp(K(k).bloque, 'ref')
        X(:, cols) = M(:, (K(k).i1 - 1)*C + (1:C));
        continue
    end
    Y = curvas.(K(k).bloque)(:, :, K(k).i1:K(k).i2);          % W x C x n
    x = reshape(log10(K(k).x(K(k).i1:K(k).i2)), 1, 1, []);
    X(:, cols) = K(k).signo * pendiente_vec(x, Y, cfg.calib_hd.min_puntos_validos);
end
end

function s = pendiente_vec(x, Y, minp)
% Pendiente de mínimos cuadrados por fila de Y (W x C x n), ignorando NaN.
ok = isfinite(Y);
x = x .* ok; Y(~ok) = 0;
n = sum(ok, 3);
sx = sum(x, 3); sy = sum(Y, 3);
sxx = sum(x.^2, 3); sxy = sum(x .* Y, 3);
s = (n .* sxy - sx .* sy) ./ (n .* sxx - sx.^2);
s(n < minp) = NaN;
end
