function C = candidatos_rango(N, cfg)
%CANDIDATOS_RANGO  Caracterizadores candidatos para la calibración.
%   C = candidatos_rango(N, cfg) enumera las pendientes de log10<R/S> vs
%   log10 m sobre un RANGO CONTIGUO de niveles diádicos m_i = floor(N/2^i),
%   para cada forma de promediar R/S. HO y HRS(p) son casos particulares:
%   con N = 512, HO (bloque >= 8) = rango 512..8 y HRS64 = rango 512..4.
%   Se exigen al menos cfg.hurst.min_escalas niveles y bloques >= 2.
%   C(k).nombre ('Re512a8': e = Ec. 13, o = rangors original), .modo,
%   .ihi, .ilo (índices de nivel, base 1), .m_max, .m_min.
ms = floor(N ./ 2.^(0:floor(log2(N))));
ms = ms(ms >= cfg.hurst.bloque_min_abs);
modos = {'ecuacion13', 'e'; 'rangors_original', 'o'};
C = struct('nombre', {}, 'modo', {}, 'ihi', {}, 'ilo', {}, 'm_max', {}, 'm_min', {});
for q = 1:size(modos, 1)
    for ihi = 1:numel(ms)
        for ilo = ihi + cfg.hurst.min_escalas - 1 : numel(ms)
            C(end+1) = struct('nombre', sprintf('R%s%da%d', modos{q,2}, ms(ihi), ms(ilo)), ...
                'modo', modos{q,1}, 'ihi', ihi, 'ilo', ilo, 'm_max', ms(ihi), 'm_min', ms(ilo)); %#ok<AGROW>
        end
    end
end
end
