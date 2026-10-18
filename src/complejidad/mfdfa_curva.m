function F = mfdfa_curva(V, ns, orden, qs)
%MFDFA_CURVA  Funciones de fluctuación multifractales F_q(n) (Kantelhardt et al., 2002).
%   F = mfdfa_curva(V, ns, orden, qs), V (C x N) -> F (C x numel(ns) x numel(qs)).
%   Perfil y = cumsum(x - media); ventanas contiguas de n muestras con
%   destendencia polinómica de grado 'orden'; varianza por ventana F2_v y
%   F_q(n) = ( media_v F2_v^(q/2) )^(1/q) (q ~= 0). h(q) = pendiente de
%   log F_q vs log n: q > 0 resalta fluctuaciones grandes, q < 0 las pequeñas;
%   el ancho h(-q) - h(q) mide la multifractalidad. Con q = 2 es el DFA.
[C, N] = size(V);
Y = cumsum(double(V) - mean(double(V), 2), 2).';
F = nan(C, numel(ns), numel(qs));
for i = 1:numel(ns)
    n = ns(i); nv = floor(N / n);
    if nv < 2 || n <= orden + 1, continue, end
    B = reshape(Y(1:nv*n, :), n, nv * C);
    t = ((0:n-1).' - (n-1)/2) / n;
    P = t .^ (0:orden);
    res = B - P * (P \ B);
    F2 = reshape(mean(res.^2, 1), nv, C);          % nv x C
    F2 = max(F2, 1e-20);
    for j = 1:numel(qs)
        F(:, i, j) = mean(F2 .^ (qs(j) / 2), 1).' .^ (1 / qs(j));
    end
end
F(F <= 1e-10) = NaN;
end
