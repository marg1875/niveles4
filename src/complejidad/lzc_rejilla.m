function S = lzc_rejilla(V, taus)
%LZC_REJILLA  Complejidad de Lempel-Ziv (LZ76) normalizada.
%   S = lzc_rejilla(V, taus), V (C x N) -> S (C x 2*numel(taus)): primero
%   binarización por la mediana y luego 3 niveles (terciles), cada una para
%   cada decimación tau (se toma una muestra de cada tau). Algoritmo de
%   Kaspar y Schuster (1987); normalización c(n) * log_a(n) / n, con a el
%   número de símbolos.
X = double(V);
C = size(X, 1);
nt = numel(taus);
S = nan(C, 2 * nt);
for it = 1:nt
    Xt = X(:, 1:taus(it):end);
    n = size(Xt, 2);
    for c = 1:C
        x = Xt(c, :);
        S(c, it) = lz76(double(x > median(x))) * log2(n) / n;
        xs = sort(x);
        q1 = xs(max(1, round(n / 3))); q2 = xs(max(1, round(2 * n / 3)));
        S(c, nt + it) = lz76((x > q1) + (x > q2)) * (log(n) / log(3)) / n;
    end
end
end

function c = lz76(s)
% Número de palabras nuevas del análisis de Lempel-Ziv 1976 (Kaspar-Schuster).
n = numel(s); c = 1; l = 1; i = 0; k = 1; kmax = 1;
while true
    if s(i + k) == s(l + k)
        k = k + 1;
        if l + k > n, c = c + 1; break, end
    else
        kmax = max(k, kmax); i = i + 1;
        if i == l
            c = c + 1; l = l + kmax;
            if l + 1 > n, break, end
            i = 0; k = 1; kmax = 1;
        else
            k = 1;
        end
    end
end
end
