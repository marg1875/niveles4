function S = fuzzyen_rejilla(V, m, rs, taus, n)
%FUZZYEN_REJILLA  Entropía difusa FuzzyEn(m, r, n, tau) (Chen et al., 2007).
%   S = fuzzyen_rejilla(V, m, rs, taus, n), V (C x N) -> S (C x numel(taus)*numel(rs)),
%   orden tau (lento), r (rápido); r relativo a std(x). Plantillas de m y
%   m + 1 puntos sin su media local; similitud exp(-(d^n)/r) con d la
%   distancia de Chebyshev; FuzzyEn = ln(phi_m) - ln(phi_m+1). Menos
%   sensible que SampEn a la elección de r (no hay umbral duro).
if nargin < 5, n = 2; end
X = double(V);
[C, N] = size(X);
S = nan(C, numel(taus) * numel(rs));
for c = 1:C
    x = X(c, :); sd = std(x);
    col = 0;
    for t = taus
        M = N - m * t;
        if M < 10 || sd < 1e-12, col = col + numel(rs); continue, end
        tri = triu(true(M), 1);
        phi = zeros(2, numel(rs));
        for q = 0:1
            mm = m + q;
            E = x((1:M).' + (0:mm-1) * t);
            E = E - mean(E, 2);
            D = zeros(M);
            for k = 1:mm, D = max(D, abs(E(:, k) - E(:, k).')); end
            dk = D(tri);
            for ir = 1:numel(rs)
                phi(q + 1, ir) = mean(exp(-(dk .^ n) / (rs(ir) * sd)));
            end
        end
        S(c, col + (1:numel(rs))) = log(phi(1, :)) - log(phi(2, :));
        col = col + numel(rs);
    end
end
end
