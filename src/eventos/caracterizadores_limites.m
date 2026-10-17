function [f, nombres] = caracterizadores_limites(V, fams, cfg)
%CARACTERIZADORES_LIMITES  Hurst por R/S en rangos de bloque definidos por límites.
%   [f, nombres] = caracterizadores_limites(V, fams, cfg), V (C x L).
%   fams: struct con nombre, modo ('ecuacion13'|'rangors_original'),
%   m_max y m_min. Se usan los niveles diádicos floor(L/2^i) que caen en
%   [m_min, m_max], de modo que el mismo caracterizador calibrado (p. ej.
%   Ro64a2 = bloques de 64 a 2, promedio de rangors) se aplica a ventanas de
%   cualquier longitud. m_max = Inf significa "desde el bloque completo".
%   f (1 x C*numel(fams)), nombres '<fam>_<canal>'.
L = size(V, 2);
ms = floor(L ./ 2.^(0:floor(log2(L))));
ms = ms(ms >= cfg.hurst.bloque_min_abs);
RS = struct('ecuacion13', rs_diadico(V, ms, 'ecuacion13'), ...
            'rangors_original', rs_diadico(V, ms, 'rangors_original'));
C = size(V, 1);
f = zeros(1, C * numel(fams));
nombres = cell(1, C * numel(fams));
for k = 1:numel(fams)
    sel = ms <= fams(k).m_max & ms >= fams(k).m_min;
    h = pendiente_hurst(ms(sel), RS.(fams(k).modo)(:, sel), cfg.hurst.min_escalas);
    f((k-1)*C + (1:C)) = h.';
    nombres((k-1)*C + (1:C)) = strcat(fams(k).nombre, '_', cfg.canales.nombres);
end
end
