function F = dfa_curva(V, ns, orden)
%DFA_CURVA  Función de fluctuación F(n) del DFA para varias escalas.
%   F = dfa_curva(V, ns, orden), V (C x N) -> F (C x numel(ns)).
%   y = cumsum(x - media) (perfil); para cada escala n con >= 1 ventana se
%   quita a cada ventana su polinomio de mínimos cuadrados de grado 'orden'
%   (1 = DFA1 lineal, como Niveles3; 2 = DFA2 cuadrático) y
%   F(n) = sqrt(sum residuos^2 / (ventanas * n)). F <= 1e-10 -> NaN.
if nargin < 3, orden = 1; end
[C, N] = size(V);
X = double(V);
Y = cumsum(X - mean(X, 2), 2).';                  % N x C
F = nan(C, numel(ns));
for i = 1:numel(ns)
    n = ns(i); nv = floor(N / n);
    if nv < 1 || n <= orden + 1, continue, end
    B = reshape(Y(1:nv*n, :), n, nv * C);
    t = ((0:n-1).' - (n-1)/2) / n;
    P = t .^ (0:orden);                            % n x (orden+1)
    res = B - P * (P \ B);                         % destendencia polinómica por ventana
    F(:, i) = sqrt(sum(reshape(res.^2, n * nv, C), 1).' / (nv * n));
end
F(F <= 1e-10) = NaN;
end
