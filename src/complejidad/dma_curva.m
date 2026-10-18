function F = dma_curva(V, ns)
%DMA_CURVA  Análisis de fluctuación con media móvil centrada (DMA; Alessio et al., 2002).
%   F = dma_curva(V, ns), V (C x N) -> F (C x numel(ns)). Perfil y; se resta
%   la media móvil centrada de n muestras (n impar) y F(n) = rms del residuo
%   lejos de los bordes. alfa = pendiente de log F vs log n (como el DFA,
%   pero sin partir en ventanas: menos sensible a la colocación de bloques).
[C, N] = size(V);
Y = cumsum(double(V) - mean(double(V), 2), 2);
F = nan(C, numel(ns));
for i = 1:numel(ns)
    n = ns(i) + mod(ns(i) + 1, 2);                 % impar
    if n >= N / 2, continue, end
    r = Y - movmean(Y, n, 2);
    h = (n - 1) / 2;
    F(:, i) = sqrt(mean(r(:, h+1:N-h).^2, 2));
end
F(F <= 1e-10) = NaN;
end
