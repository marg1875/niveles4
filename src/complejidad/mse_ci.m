function S = mse_ci(V, Ss, m, r)
%MSE_CI  Índice de complejidad de la entropía multiescala (Costa et al., 2002).
%   S = mse_ci(V, Ss, m, r), V (C x N) -> S (C x numel(Ss)): suma de
%   SampEn(m, r) de las series de grano grueso (medias de bloques de tamaño
%   e = 1..S) para cada S en Ss. r es relativo a la std de la serie
%   ORIGINAL, como en Costa et al. Las escalas con < 10 plantillas se omiten.
X = double(V);
[C, N] = size(X);
Smax = max(Ss);
sd = std(X, 0, 2);
se = nan(C, Smax);
for e = 1:Smax
    nb = floor(N / e);
    if nb - m < 10, break, end
    G = reshape(mean(reshape(X(:, 1:nb*e).', e, nb * C), 1), nb, C).';
    for c = 1:C
        rel = r * sd(c) / max(std(G(c, :)), eps);   % misma r absoluta que la serie original
        se(c, e) = sampen_rejilla(G(c, :), m, rel, 1);
    end
end
S = nan(C, numel(Ss));
for i = 1:numel(Ss), S(:, i) = sum(se(:, 1:Ss(i)), 2, 'omitnan'); end
end
