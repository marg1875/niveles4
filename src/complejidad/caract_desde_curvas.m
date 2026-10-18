function [X, nombres] = caract_desde_curvas(M, K, lay, cfg)
%CARACT_DESDE_CURVAS  Candidatos de caract_registro a partir de las filas de caract_curvas.
%   [X, nombres] = caract_desde_curvas(M, K, lay, cfg)
%   M (ventanas x lay.ancho). X (ventanas x numel(K)*C), nombres '<cand>_<canal>'.
%     'pend' pendiente de mínimos cuadrados de log10(curva) vs log10(x) en
%            i1:i2 (puntos NaN fuera; NaN si quedan < 3), por signo
%     'dif'  pendiente(bloque) - pendiente(bloque2) (ancho multifractal)
%     'esc'  valor directo (columna i1 del bloque)
C = lay.C; W = size(M, 1);
X = nan(W, numel(K) * C);
nombres = cell(1, numel(K) * C);
idx = containers.Map({lay.bloques.nombre}, num2cell(1:numel(lay.bloques)));
cache = struct();
for k = 1:numel(K)
    cols = (k-1)*C + (1:C);
    nombres(cols) = strcat(K(k).nombre, '_', cfg.canales.nombres);
    b = lay.bloques(idx(K(k).bloque));
    if strcmp(K(k).tipo, 'esc')
        X(:, cols) = M(:, b.o + (K(k).i1 - 1)*C + (1:C));
        continue
    end
    s = pendiente(curva(b), K(k));
    if strcmp(K(k).tipo, 'dif')
        s = s - pendiente(curva(lay.bloques(idx(K(k).bloque2))), K(k));
    end
    X(:, cols) = K(k).signo * s;
end

    function Y = curva(b)
        if ~isfield(cache, b.nombre)
            cache.(b.nombre) = log10(reshape(M(:, b.o + (1:b.K*C)), [], C, b.K));
        end
        Y = cache.(b.nombre);
    end
end

function s = pendiente(Y, k)
Y = Y(:, :, k.i1:k.i2);
x = reshape(log10(k.x(k.i1:k.i2)), 1, 1, []);
ok = isfinite(Y);
x = x .* ok; Y(~ok) = 0;
n = sum(ok, 3);
sx = sum(x, 3); sy = sum(Y, 3);
s = (n .* sum(x .* Y, 3) - sx .* sy) ./ (n .* sum(x.^2, 3) - sx.^2);
s(n < 3) = NaN;
end
