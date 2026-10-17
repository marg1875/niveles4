function [H, nvalidos] = pendiente_hurst(ms, RS, min_escalas)
%PENDIENTE_HURST  Exponente de Hurst = pendiente de log10<R/S> vs log10 m.
%   [H, nvalidos] = pendiente_hurst(ms, RS, min_escalas), con RS (C x numel(ms)).
%   Ecs. 14-15 del artículo: <R/S>_m ~ c m^H  ->  log10<R/S> = H log10 m + b.
%   Si una fila tiene menos de min_escalas puntos válidos, H = NaN.
C = size(RS, 1);
H = nan(C, 1);
nvalidos = zeros(C, 1);
x = log10(ms(:));
for c = 1:C
    ok = isfinite(RS(c,:)).' & RS(c,:).' > 0;
    nvalidos(c) = nnz(ok);
    if nvalidos(c) >= min_escalas
        xc = x(ok); yc = log10(RS(c,ok)).';
        xm = mean(xc);
        H(c) = sum((xc - xm) .* (yc - mean(yc))) / sum((xc - xm).^2);
    end
end
end
