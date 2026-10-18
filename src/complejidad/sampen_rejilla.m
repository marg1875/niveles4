function S = sampen_rejilla(V, ms, rs, taus)
%SAMPEN_REJILLA  Entropía muestral SampEn(m, r, tau) para una rejilla de parámetros.
%   S = sampen_rejilla(V, ms, rs, taus), V (C x N) -> S (C x numel(taus)*numel(ms)*numel(rs))
%   en el orden: tau (más lento), m, r (más rápido). r es relativo a std(x).
%   Richman y Moorman (2000): B = pares de plantillas de longitud m (paso
%   tau) a distancia de Chebyshev <= r, A = ídem con m + 1, SampEn = -ln(A/B).
%   Se usan las mismas N - m_max*tau plantillas para todo m (comparables).
%   Las distancias se acumulan dimensión a dimensión y se cuentan para
%   todos los r a la vez. A = 0 -> NaN (indefinida).
X = double(V);
[C, N] = size(X);
mmax = max(ms);
S = nan(C, numel(taus) * numel(ms) * numel(rs));
for c = 1:C
    x = X(c, :); sd = std(x);
    col = 0;
    for t = taus
        M = N - mmax * t;
        if M < 10 || sd < 1e-12, col = col + numel(ms) * numel(rs); continue, end
        cnt = zeros(mmax + 1, numel(rs));
        D = zeros(M);
        tri = triu(true(M), 1);
        for k = 0:mmax
            v = x((1:M) + k * t);
            D = max(D, abs(v.' - v));
            dk = D(tri);
            for ir = 1:numel(rs), cnt(k + 1, ir) = nnz(dk <= rs(ir) * sd); end
        end
        for m = ms
            for ir = 1:numel(rs)
                col = col + 1;
                B = cnt(m, ir); A = cnt(m + 1, ir);
                if A > 0 && B > 0, S(c, col) = -log(A / B); end
            end
        end
    end
end
end
