function B = boxcount_curva(V, ss)
%BOXCOUNT_CURVA  Conteo de cajas de la gráfica de la señal (dimensión box-counting).
%   B = boxcount_curva(V, ss), V (C x N) -> B (C x numel(ss)). La gráfica se
%   lleva al cuadrado unitario (tiempo y amplitud normalizados) y se cubre
%   con una rejilla de s x s cajas: en cada columna se cuentan las cajas que
%   toca la curva entre su mínimo y su máximo. D = pendiente de log B vs log s.
[C, N] = size(V);
X = double(V);
X = (X - min(X, [], 2)) ./ max(max(X, [], 2) - min(X, [], 2), eps);
B = nan(C, numel(ss));
for i = 1:numel(ss)
    s = ss(i);
    if s > N - 1, continue, end
    b = floor((0:N-1) / (N - 1) * s); b(end) = s - 1;   % columna de cada muestra
    tot = zeros(C, 1);
    for k = 0:s-1
        idx = find(b == k);
        if isempty(idx), continue, end
        idx = [max(idx(1) - 1, 1), idx];                % une con la columna anterior
        lo = floor(min(X(:, idx), [], 2) * s); hi = floor(min(max(X(:, idx), [], 2), 1 - eps) * s);
        tot = tot + (hi - lo + 1);
    end
    B(:, i) = tot;
end
end
