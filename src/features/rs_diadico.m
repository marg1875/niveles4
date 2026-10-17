function RS = rs_diadico(V, ms, modo)
%RS_DIADICO  Rango reescalado medio <R/S>_m para varios tamaños de bloque.
%   RS = rs_diadico(V, ms), con V (C x N: una serie por fila) y ms vector de
%   tamaños de bloque. RS es C x numel(ms).
%   RS = rs_diadico(V, ms, 'rangors_original') reproduce exactamente
%   rangors.m, que escribe mean(R / S) con '/' MATRICIAL: el resultado es
%   el cociente de mínimos cuadrados sum(R.*S)/sum(S.^2), no la media de
%   R./S de la Ec. 13. Por defecto ('ecuacion13') se usa la Ec. 13.
%
%   Port vectorizado del bucle interno de rangors.m (grupo Benavides-Bravo /
%   Martinez-Peon; Ecs. 10-13 del artículo). Para cada tamaño m, la serie se
%   parte en d = floor(N/m) subgrupos contiguos L_i (se descarta el resto);
%   en cada uno:  G = L_i - media(L_i);  V = cumsum(G);  R = max V - min V;
%   S = std(L_i) (normalización N-1, como el original). <R/S>_m = media de
%   R/S sobre los subgrupos. Diferencia documentada: si S = 0 el subgrupo se
%   omite (el original devolvería NaN/Inf).
if nargin < 3, modo = 'ecuacion13'; end
original = strcmp(modo, 'rangors_original');
[C, N] = size(V);
RS = nan(C, numel(ms));
Vt = double(V.');                                  % N x C
for i = 1:numel(ms)
    m = ms(i);
    nb = floor(N / m);
    if m < 2 || nb < 1, continue, end
    B = reshape(Vt(1:nb*m, :), m, nb, C);           % m x nb x C
    Y = cumsum(B - mean(B, 1), 1);
    R = max(Y, [], 1) - min(Y, [], 1);
    S = std(B, 0, 1);
    if original
        RS(:, i) = reshape(sum(R.*S, 2) ./ sum(S.^2, 2), C, 1);
    else
        q = R ./ S;
        q(S < 1e-12) = NaN;
        RS(:, i) = reshape(mean(q, 2, 'omitnan'), C, 1);
    end
end
end
