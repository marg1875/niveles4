function [Gv, Ga] = agregados_curva(V, ms)
%AGREGADOS_CURVA  Varianza agregada y momento absoluto de las medias por bloque.
%   [Gv, Ga] = agregados_curva(V, ms), V (C x N) -> (C x numel(ms)) cada una.
%   Para cada tamaño m se promedia x en bloques contiguos: Gv(m) = varianza
%   de esas medias y Ga(m) = media de |media_bloque - media|. Para un
%   proceso con exponente H: Gv ~ m^(2H-2) y Ga ~ m^(H-1) (Taqqu et al., 1995).
[C, N] = size(V);
X = double(V).';
Gv = nan(C, numel(ms)); Ga = Gv;
for i = 1:numel(ms)
    m = ms(i); nb = floor(N / m);
    if nb < 4, continue, end
    M = squeeze(mean(reshape(X(1:nb*m, :), m, nb, C), 1));   % nb x C
    if C == 1, M = M(:); end
    Gv(:, i) = var(M, 1, 1).';
    Ga(:, i) = mean(abs(M - mean(M, 1)), 1).';
end
Gv(Gv <= 1e-20) = NaN; Ga(Ga <= 1e-20) = NaN;
end
