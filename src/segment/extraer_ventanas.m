function V = extraer_ventanas(X, centros, L)
%EXTRAER_VENTANAS  Ventanas de L muestras centradas en cada ancla.
%   V = extraer_ventanas(X, centros, L) con X (C x N). Devuelve V (C x L x E).
%   La ventana de la ancla c abarca [c - L/2, c + L/2 - 1] (L par) o
%   [c - floor(L/2), c + floor(L/2)] (L impar); con L = 512 coincide
%   exactamente con la ranura de la ancla. Si la ventana se saldría del
%   registro (L grande cerca de los bordes) se desplaza hasta quedar dentro.
E = numel(centros);
V = zeros(size(X,1), L, E, 'like', X);
for e = 1:E
    s = min(max(centros(e) - floor(L/2), 1), size(X, 2) - L + 1);   % se desplaza dentro del registro si se sale
    V(:,:,e) = X(:, s:s+L-1);
end
end
