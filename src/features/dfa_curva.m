function F = dfa_curva(V, ns, orden, paso)
%DFA_CURVA  Función de fluctuación F(n) del DFA para varias escalas.
%   F = dfa_curva(V, ns, orden, paso), V (C x N) -> F (C x numel(ns)).
%   y = cumsum(x - media) (perfil); para cada escala n se quita a cada
%   ventana su polinomio de mínimos cuadrados de grado 'orden' (1 = DFA1
%   lineal, como Niveles3; 2 = DFA2; 3 = DFA3) y
%   F(n) = sqrt(media de los residuos^2).
%   paso (opcional, fracción de n): 1 = ventanas contiguas sin traslape
%   (por defecto, Niveles3); 0.5 = ventanas deslizantes con 50 % de
%   traslape (más ventanas por escala: menos varianza en escalas grandes).
%   Escalas con menos de 1 ventana o n <= orden + 1 -> NaN; F <= 1e-10 -> NaN.
if nargin < 3, orden = 1; end
if nargin < 4, paso = 1; end
[C, N] = size(V);
X = double(V);
Y = cumsum(X - mean(X, 2), 2).';                  % N x C
F = nan(C, numel(ns));
for i = 1:numel(ns)
    n = ns(i);
    if n > N || n <= orden + 1, continue, end
    d = max(1, round(paso * n));
    ini = 1:d:N-n+1;
    if paso >= 1, ini = 1:n:floor(N/n)*n - n + 1; end
    nv = numel(ini);
    idx = ini + (0:n-1).';                         % n x nv
    B = reshape(Y(idx(:), :), n, nv * C);
    t = ((0:n-1).' - (n-1)/2) / n;
    P = t .^ (0:orden);
    res = B - P * (P \ B);                         % destendencia polinómica por ventana
    F(:, i) = sqrt(sum(reshape(res.^2, n * nv, C), 1).' / (nv * n));
end
F(F <= 1e-10) = NaN;
end
