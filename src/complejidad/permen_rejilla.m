function S = permen_rejilla(V, ds, taus)
%PERMEN_REJILLA  Entropía de permutación (Bandt-Pompe) y su versión ponderada.
%   S = permen_rejilla(V, ds, taus), V (C x N) -> S (C x 2*numel(ds)*numel(taus)):
%   primero PE y luego WPE (Fadlallah et al., 2013), cada una en el orden
%   d (lento), tau (rápido). Entropía de Shannon de los patrones ordinales
%   de d puntos con retardo tau, normalizada por log(d!). WPE pondera cada
%   patrón por la varianza de su vector (resta peso a los tramos planos).
%   Combinaciones con menos de d! plantillas -> NaN (estimación sesgada).
X = double(V);
[C, N] = size(X);
nd = numel(ds); nt = numel(taus);
S = nan(C, 2 * nd * nt);
for id = 1:nd
    d = ds(id);
    pot = d .^ (0:d-1);
    for it = 1:nt
        t = taus(it);
        M = N - (d - 1) * t;
        if M < factorial(d), continue, end
        k = (id - 1) * nt + it;
        for c = 1:C
            x = X(c, :); E = x((1:M).' + (0:d-1) * t);
            [~, o] = sort(E, 2);
            [~, ~, g] = unique((o - 1) * pot.');
            p = accumarray(g, 1) / M;
            w = var(E, 1, 2);
            pw = accumarray(g, w) / max(sum(w), eps);
            S(c, k) = -sum(p .* log(p)) / log(factorial(d));
            pw = pw(pw > 0);
            S(c, nd * nt + k) = -sum(pw .* log(pw)) / log(factorial(d));
        end
    end
end
end
