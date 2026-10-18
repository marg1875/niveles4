function E = wavelet_curva(V, wname, J)
%WAVELET_CURVA  Energía media de los detalles wavelet por octava (Abry-Veitch).
%   E = wavelet_curva(V, wname, J), V (C x N) -> E (C x J): E(j) = media de
%   d_j^2 (coeficientes de detalle del nivel j, wavedec). log2 E(j) ~
%   (2H + 1) j: la pendiente sobre un rango de octavas estima H.
%   Niveles no calculables (N demasiado corto para el filtro) -> NaN.
[C, N] = size(V);
E = nan(C, J);
Jmax = min(J, wmaxlev(N, wname));
if Jmax < 1, return, end
for c = 1:C
    [cf, l] = wavedec(double(V(c, :)), Jmax, wname);
    for j = 1:Jmax
        d = detcoef(cf, l, j);
        E(c, j) = mean(d.^2);
    end
end
E(E <= 1e-20) = NaN;
end
