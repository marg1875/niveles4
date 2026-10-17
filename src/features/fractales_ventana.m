function [f, nombres] = fractales_ventana(V, cfg, base)
%FRACTALES_VENTANA  Los 8 algoritmos fractales de Niveles3 + el par calibrado.
%   [f, nombres] = fractales_ventana(V, cfg, base), V (C x L).
%   Familias (cfg.fractales.familias), cada una con un valor por canal:
%     RS    Hurst R/S con destendencia (hurst_rs_tendencia)
%     HIG   dimensión de Higuchi (higuchi_dimension)
%     DFA   exponente alfa (dfa_exponente)
%     VAR   pendiente/2 del semivariograma (variograma_pendiente)
%     PROM  media de RS, HIG, DFA y VAR (como Niveles3; mezcla D con H)
%     HO    Hurst original diádico (hurst_ho)
%     HRS64 Hurst con 64 particiones (hurst_hrs)
%     HV    R/S del semivariograma (hurst_hv, script del grupo)
%   base (opcional): struct de familias de caracterizadores_limites (p. ej.
%   el par calibrado Re16a2 + Re64a4); se colocan al inicio.
%   f (1 x C*nfam), nombres '<familia>_<canal>'.
if nargin < 3, base = struct('nombre', {}); end
ch = cfg.canales.nombres;
fb = []; nb = {};
if ~isempty(base), [fb, nb] = caracterizadores_limites(V, base, cfg); end
rs  = hurst_rs_tendencia(V, cfg);
hig = higuchi_dimension(V, cfg);
dfa = dfa_exponente(V, cfg);
vr  = variograma_pendiente(V, cfg);
M = struct('RS', rs, 'HIG', hig, 'DFA', dfa, 'VAR', vr, ...
    'PROM', mean([rs hig dfa vr], 2, 'omitnan'), ...
    'HO', hurst_ho(V, cfg), 'HRS64', hurst_hrs(V, cfg.hurst.p_hrs_principal, cfg), ...
    'HV', hurst_hv(V, cfg));
fam = cfg.fractales.familias;
f = zeros(1, numel(ch) * numel(fam));
nombres = cell(1, numel(f));
for k = 1:numel(fam)
    f((k-1)*numel(ch) + (1:numel(ch))) = M.(fam{k}).';
    nombres((k-1)*numel(ch) + (1:numel(ch))) = strcat(fam{k}, '_', ch);
end
f = [fb, f];
nombres = [nb, nombres];
end
