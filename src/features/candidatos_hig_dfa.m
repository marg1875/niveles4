function [K, lay] = candidatos_hig_dfa(cfg, baseH, L)
%CANDIDATOS_HIG_DFA  Rejilla de calibración de Higuchi y DFA (etapas 10 y 12).
%   [K, lay] = candidatos_hig_dfa(cfg, baseH, L)
%   Higuchi: D = -pendiente de log L(k) en k = kmin..kmax (todos los k
%     enteros), kmin en cfg.calib_hd.hig_kmin, kmax en cfg.calib_hd.hig_kmax
%     (y, si se da L, kmax <= floor(L/4), como Niveles3), con al menos
%     cfg.calib_hd.hig_min_puntos valores de k. Nombre 'HIG<kmin>a<kmax>'.
%   DFA: alfa = pendiente de log F(n) sobre un rango contiguo de la rejilla
%     cfg.calib_hd.dfa_ns (>= cfg.calib_hd.dfa_min_escalas), para cada
%     orden de destendencia (cfg.calib_hd.dfa_ordenes: 1 'l', 2 'c', 3 't')
%     y paso de ventana (cfg.calib_hd.dfa_pasos: 1 sin traslape, 0.5 con
%     50 % 'o'). Con L: n <= floor(L/3) sin traslape y n <= floor(L/2) con
%     traslape. Nombre 'DFA<orden>[o]<nmin>a<nmax>', p. ej. 'DFAc4a24', 'DFAto8a64'.
%   Referencias: 'HIGn3', 'DFAn3' (Niveles3) y las familias de Hurst de baseH.
%   cfg.calib_hd.prefijo (opcional, p. ej. 'W250') se antepone a todos los nombres.
%   K: struct con nombre, tipo ('HIG'|'DFA'|'REF'), bloque, i1:i2, signo y x.
%   lay: disposición de la fila de curvas_ventana (lay.bloques).
p = cfg.calib_hd;
if nargin < 3, L = Inf; end
ordenes = 1:2; if isfield(p, 'dfa_ordenes'), ordenes = p.dfa_ordenes; end
pasos = 1;     if isfield(p, 'dfa_pasos'), pasos = p.dfa_pasos; end
pre = '';      if isfield(p, 'prefijo'), pre = p.prefijo; end
C = numel(cfg.canales.nombres);
kmax = p.hig_kmax(p.hig_kmax <= floor(L / 4));
lay.C = C; lay.kmax = max(kmax); lay.ns = p.dfa_ns(p.dfa_ns <= floor(L / 2));
lay.ordenes = ordenes; lay.pasos = pasos;
lay.nref = numel(baseH) + 2;
o = lay.nref * C;
lay.bloques = struct('nombre', 'hig', 'o', o, 'K', lay.kmax);
o = o + lay.kmax * C;
for a = ordenes
    for b = 1:numel(pasos)
        lay.bloques(end+1) = struct('nombre', sprintf('dfa%d_%d', a, b), 'o', o, 'K', numel(lay.ns));
        o = o + numel(lay.ns) * C;
    end
end
lay.ancho = o;
K = struct('nombre', {}, 'tipo', {}, 'bloque', {}, 'i1', {}, 'i2', {}, 'signo', {}, 'x', {});
refs = [{baseH.nombre}, {'HIGn3', 'DFAn3'}];
for r = 1:numel(refs)
    K(end+1) = struct('nombre', [pre refs{r}], 'tipo', 'REF', 'bloque', 'ref', 'i1', r, 'i2', r, 'signo', 1, 'x', []); %#ok<AGROW>
end
for a = p.hig_kmin
    for b = kmax
        if b - a + 1 < p.hig_min_puntos, continue, end
        K(end+1) = struct('nombre', sprintf('%sHIG%da%d', pre, a, b), 'tipo', 'HIG', 'bloque', 'hig', ...
            'i1', a, 'i2', b, 'signo', -1, 'x', 1:lay.kmax); %#ok<AGROW>
    end
end
ns = lay.ns;
letra = 'lct';
for a = ordenes
    for b = 1:numel(pasos)
        nmax = floor(L / 3); suf = '';
        if pasos(b) < 1, nmax = floor(L / 2); suf = 'o'; end
        for i = 1:numel(ns)
            if ns(i) <= a + 1, continue, end
            for j = i + p.dfa_min_escalas - 1 : numel(ns)
                if ns(j) > nmax, break, end
                K(end+1) = struct('nombre', sprintf('%sDFA%s%s%da%d', pre, letra(a), suf, ns(i), ns(j)), ...
                    'tipo', 'DFA', 'bloque', sprintf('dfa%d_%d', a, b), 'i1', i, 'i2', j, 'signo', 1, 'x', ns); %#ok<AGROW>
            end
        end
    end
end
end
