function [K, lay] = caract_registro(cfg, L, baseH)
%CARACT_REGISTRO  Registro de familias de caracterización y sus rejillas (etapa 13).
%   [K, lay] = caract_registro(cfg, L, baseH) para ventanas de L muestras.
%   lay.bloques: qué calcula caract_curvas por ventana (curvas o valores).
%   K: un candidato por fila con
%     nombre   'W<L><familia><parámetros>' (sin '_': el canal va después)
%     familia  grupo para los rankings (HIG, DFA1, MFh+3, PE, SAMP, ...)
%     tipo     'pend' (pendiente log-log en un rango), 'dif' (diferencia de
%              dos pendientes, ancho multifractal) o 'esc' (valor directo)
%     bloque, i1, i2 (y bloque2 para 'dif'), x (abscisas), signo
%   Rejillas en cfg.complejidad (ver config_niveles.m). Los rangos de
%   escalas se toman de una rejilla dispersa (nmin x nmax) para que todas
%   las familias tengan un número de candidatos comparable.
p = cfg.complejidad;
C = numel(cfg.canales.nombres);
pre = sprintf('W%d', L);
lay = struct('C', C, 'L', L, 'bloques', struct('nombre', {}, 'tipo', {}, 'K', {}, 'o', {}, 'arg', {}));
K = struct('nombre', {}, 'familia', {}, 'tipo', {}, 'bloque', {}, 'bloque2', {}, 'i1', {}, 'i2', {}, 'x', {}, 'signo', {});

% --- referencias (Hurst calibrado de etapas previas, Higuchi/DFA de Niveles3)
nref = numel(baseH) + 2;
lay = agregar(lay, 'ref', 'esc', nref, struct('baseH', baseH));
nomRef = [{baseH.nombre}, {'HIGn3', 'DFAn3'}];
for r = 1:nref, K = esc(K, [pre nomRef{r}], 'REF', 'ref', r); end

% --- Higuchi: k = kmin..kmax (todos los k)
kmax = p.hig_kmax(p.hig_kmax <= floor(L / 4));
lay = agregar(lay, 'hig', 'curva', max(kmax), struct());
for a = p.hig_kmin
    for b = kmax
        if b - a + 1 >= 5
            K = pend(K, sprintf('%sHIG%da%d', pre, a, b), 'HIG', 'hig', a, b, 1:max(kmax), -1);
        end
    end
end

% --- DFA (orden 1-3, con y sin traslape), DMA y MFDFA sobre la rejilla de escalas
ns = p.ns(p.ns <= floor(L / 3));
R = rangos(ns, p.n_min, p.n_max, 4);
letra = 'lct';
for o = p.dfa_ordenes
    for ip = 1:numel(p.dfa_pasos)
        b = sprintf('dfa%d_%d', o, ip); suf = ''; if p.dfa_pasos(ip) < 1, suf = 'o'; end
        lay = agregar(lay, b, 'curva', numel(ns), struct('ns', ns, 'orden', o, 'paso', p.dfa_pasos(ip)));
        for r = R.', if ns(r(1)) > o + 1
            K = pend(K, sprintf('%sDFA%s%s%da%d', pre, letra(o), suf, ns(r(1)), ns(r(2))), ['DFA' num2str(o)], b, r(1), r(2), ns, 1);
        end, end
    end
end
nsd = ns(ns >= 5);
lay = agregar(lay, 'dma', 'curva', numel(nsd), struct('ns', nsd));
for r = rangos(nsd, p.n_min, p.n_max, 4).'
    K = pend(K, sprintf('%sDMA%da%d', pre, nsd(r(1)), nsd(r(2))), 'DMA', 'dma', r(1), r(2), nsd, 1);
end
for o = p.mf_ordenes
    for iq = 1:numel(p.mf_q)
        b = sprintf('mf%d_%d', o, iq);
        lay = agregar(lay, b, 'curva', numel(ns), struct('ns', ns, 'orden', o, 'iq', iq));
    end
    for r = R.'
        if ns(r(1)) <= o + 1, continue, end
        tag = sprintf('o%d_%da%d', o, ns(r(1)), ns(r(2)));
        for iq = 1:numel(p.mf_q)
            q = p.mf_q(iq);
            K = pend(K, strrep(sprintf('%sMFh%s%d%s', pre, signo_txt(q), abs(q), tag), '_', 'n'), ...
                sprintf('MFh%s%d', signo_txt(q), abs(q)), sprintf('mf%d_%d', o, iq), r(1), r(2), ns, 1);
        end
        for w = p.mf_anchos
            i1 = find(p.mf_q == -w); i2 = find(p.mf_q == w);
            if isempty(i1) || isempty(i2), continue, end
            K(end+1) = struct('nombre', strrep(sprintf('%sMFw%d%s', pre, w, tag), '_', 'n'), 'familia', 'MFancho', ...
                'tipo', 'dif', 'bloque', sprintf('mf%d_%d', o, i1), 'bloque2', sprintf('mf%d_%d', o, i2), ...
                'i1', r(1), 'i2', r(2), 'x', ns, 'signo', 1); %#ok<AGROW>
        end
    end
end

% --- varianza agregada y momentos absolutos
ms = p.ms(p.ms <= floor(L / 4));
lay = agregar(lay, 'agv', 'curva', numel(ms), struct('ms', ms));
lay = agregar(lay, 'abm', 'curva', numel(ms), struct('ms', ms));
for r = rangos(ms, ms, ms, 4).'
    K = pend(K, sprintf('%sAGV%da%d', pre, ms(r(1)), ms(r(2))), 'AGV', 'agv', r(1), r(2), ms, 1);
    K = pend(K, sprintf('%sABM%da%d', pre, ms(r(1)), ms(r(2))), 'ABM', 'abm', r(1), r(2), ms, 1);
end

% --- Hurst por wavelets (octavas j1..j2)
J = floor(log2(L)) - 2;
for iw = 1:numel(p.wavelets)
    w = p.wavelets{iw}; b = ['wav' num2str(iw)];
    lay = agregar(lay, b, 'curva', J, struct('wname', w, 'J', J));
    for j1 = 1:J
        for j2 = j1 + 2:J
            K = pend(K, sprintf('%sWV%sj%da%d', pre, strrep(w, 'haar', 'hr'), j1, j2), 'WAV', b, j1, j2, 2 .^ (1:J), 1);
        end
    end
end

% --- exponente aperiódico 1/f (Welch)
nfft = min(256, L);
f = (0:floor(nfft/2)) * cfg.fs / nfft;
lay = agregar(lay, 'psd', 'curva', numel(f), struct('nfft', nfft, 'f', f));
for a = p.psd_f
    for b = p.psd_f
        if b < 2 * a, continue, end
        i1 = find(f >= a, 1); i2 = find(f <= b, 1, 'last');
        if isempty(i1) || i2 - i1 + 1 < 4, continue, end
        K = pend(K, sprintf('%sPSD%da%d', pre, a, b), 'PSD', 'psd', i1, i2, max(f, eps), -1);
    end
end

% --- box-counting
ss = p.box_s(p.box_s <= L / 2);
lay = agregar(lay, 'box', 'curva', numel(ss), struct('ss', ss));
for r = rangos(ss, ss, ss, 3).'
    K = pend(K, sprintf('%sBOX%da%d', pre, ss(r(1)), ss(r(2))), 'BOX', 'box', r(1), r(2), ss, 1);
end

% --- dimensiones directas y Hjorth
lay = agregar(lay, 'dir', 'esc', 7, struct());
nd = {'KATZ', 'PETd', 'PETm', 'SEV', 'HJmob', 'HJcom', 'LL'};
fd = {'KATZ', 'PET', 'PET', 'SEV', 'HJORTH', 'HJORTH', 'LL'};
for i = 1:7, K = esc(K, [pre nd{i}], fd{i}, 'dir', i); end

% --- entropías y complejidad
lay = agregar(lay, 'samp', 'esc', numel(p.se_tau) * numel(p.se_m) * numel(p.se_r), struct());
i = 0;
for t = p.se_tau, for m = p.se_m, for r = p.se_r
    i = i + 1; K = esc(K, sprintf('%sSEm%dr%dt%d', pre, m, round(100 * r), t), 'SAMPEN', 'samp', i);
end, end, end
lay = agregar(lay, 'fuz', 'esc', numel(p.fe_tau) * numel(p.fe_r), struct());
i = 0;
for t = p.fe_tau, for r = p.fe_r
    i = i + 1; K = esc(K, sprintf('%sFEr%dt%d', pre, round(100 * r), t), 'FUZZYEN', 'fuz', i);
end, end
lay = agregar(lay, 'pe', 'esc', 2 * numel(p.pe_d) * numel(p.pe_tau), struct());
i = 0;
for v = {'PE', 'WPE'}, for d = p.pe_d, for t = p.pe_tau
    i = i + 1; K = esc(K, sprintf('%s%sd%dt%d', pre, v{1}, d, t), v{1}, 'pe', i);
end, end, end
lay = agregar(lay, 'lz', 'esc', 2 * numel(p.lz_tau), struct());
i = 0;
for v = {'LZ2', 'LZ3'}, for t = p.lz_tau
    i = i + 1; K = esc(K, sprintf('%s%st%d', pre, v{1}, t), 'LZC', 'lz', i);
end, end
lay = agregar(lay, 'mse', 'esc', numel(p.mse_S), struct());
for i = 1:numel(p.mse_S), K = esc(K, sprintf('%sMSE%d', pre, p.mse_S(i)), 'MSE', 'mse', i); end
end

%% ----------------------------------------------------------------------
function lay = agregar(lay, nombre, tipo, nK, arg)
o = 0;
if ~isempty(lay.bloques), o = lay.bloques(end).o + lay.bloques(end).K * lay.C; end
lay.bloques(end+1) = struct('nombre', nombre, 'tipo', tipo, 'K', nK, 'o', o, 'arg', arg);
lay.ancho = o + nK * lay.C;
end

function K = pend(K, nombre, fam, bloque, i1, i2, x, signo)
K(end+1) = struct('nombre', nombre, 'familia', fam, 'tipo', 'pend', 'bloque', bloque, 'bloque2', '', ...
    'i1', i1, 'i2', i2, 'x', x, 'signo', signo);
end

function K = esc(K, nombre, fam, bloque, i)
K(end+1) = struct('nombre', nombre, 'familia', fam, 'tipo', 'esc', 'bloque', bloque, 'bloque2', '', ...
    'i1', i, 'i2', i, 'x', [], 'signo', 1);
end

function R = rangos(x, nmins, nmaxs, minp)
% Pares de índices (i1, i2) con x(i1) en nmins, x(i2) en nmaxs y >= minp puntos.
R = zeros(0, 2);
for a = find(ismember(x, nmins))
    for b = find(ismember(x, nmaxs))
        if b - a + 1 >= minp, R(end+1, :) = [a b]; end %#ok<AGROW>
    end
end
end

function s = signo_txt(q)
if q < 0, s = 'm'; else, s = 'p'; end
end
