function [K, lay] = candidatos_hig_dfa(cfg, baseH)
%CANDIDATOS_HIG_DFA  Rejilla de calibración de Higuchi y DFA (etapa 10).
%   [K, lay] = candidatos_hig_dfa(cfg, baseH)
%   Higuchi: D = -pendiente de log L(k) en k = kmin..kmax (todos los k
%     enteros), kmin en cfg.calib_hd.hig_kmin, kmax en cfg.calib_hd.hig_kmax,
%     al menos cfg.calib_hd.hig_min_puntos valores de k. Nombre 'HIG<kmin>a<kmax>'.
%   DFA: alfa = pendiente de log F(n) sobre un rango contiguo de la rejilla
%     de escalas cfg.calib_hd.dfa_ns (>= cfg.calib_hd.dfa_min_escalas), con
%     destendencia lineal (DFA1, 'DFAl<nmin>a<nmax>') o cuadrática (DFA2, 'DFAc...').
%   Referencias: 'HIGn3' y 'DFAn3' (parámetros por defecto de Niveles3) y
%     las familias de Hurst de baseH (caracterizadores_limites).
%   K: struct con nombre, tipo ('HIG'|'DFA'|'REF'), bloque (curva o
%   columna de referencia), rango de índices i1:i2, signo y x (abscisas).
%   lay: disposición de la fila de curvas_ventana.
p = cfg.calib_hd;
C = numel(cfg.canales.nombres);
lay.C = C; lay.kmax = max(p.hig_kmax); lay.ns = p.dfa_ns;
lay.nref = numel(baseH) + 2;                                   % familias de referencia
lay.o_ref = 0;
lay.o_hig = lay.nref * C;
lay.o_dfa1 = lay.o_hig + lay.kmax * C;
lay.o_dfa2 = lay.o_dfa1 + numel(lay.ns) * C;
lay.ancho = lay.o_dfa2 + numel(lay.ns) * C;
K = struct('nombre', {}, 'tipo', {}, 'bloque', {}, 'i1', {}, 'i2', {}, 'signo', {}, 'x', {});
refs = [{baseH.nombre}, {'HIGn3', 'DFAn3'}];
for r = 1:numel(refs)
    K(end+1) = struct('nombre', refs{r}, 'tipo', 'REF', 'bloque', 'ref', 'i1', r, 'i2', r, 'signo', 1, 'x', []); %#ok<AGROW>
end
for a = p.hig_kmin
    for b = p.hig_kmax
        if b - a + 1 < p.hig_min_puntos, continue, end
        K(end+1) = struct('nombre', sprintf('HIG%da%d', a, b), 'tipo', 'HIG', 'bloque', 'hig', ...
            'i1', a, 'i2', b, 'signo', -1, 'x', 1:lay.kmax); %#ok<AGROW>
    end
end
ns = lay.ns;
for o = 1:2
    letra = 'lc';
    for i = 1:numel(ns)
        if ns(i) <= o + 1, continue, end
        for j = i + p.dfa_min_escalas - 1 : numel(ns)
            K(end+1) = struct('nombre', sprintf('DFA%s%da%d', letra(o), ns(i), ns(j)), 'tipo', 'DFA', ...
                'bloque', sprintf('dfa%d', o), 'i1', i, 'i2', j, 'signo', 1, 'x', ns); %#ok<AGROW>
        end
    end
end
end
