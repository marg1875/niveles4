function [f, nesc] = extraer_hurst_ventana(V, cfg)
%EXTRAER_HURST_VENTANA  HO, HRS(p) y HV de todos los canales de una ventana.
%   [f, nesc] = extraer_hurst_ventana(V, cfg), V (C x L). Devuelve f
%   (1 x C*nmet) en el orden de nombres_caracteristicas(): primero HO de
%   todos los canales, luego HRS<p> para cada p, y al final HV. nesc tiene
%   el mismo orden e indica cuántas escalas entraron en cada regresión.
%
%   Los niveles diádicos de R/S se calculan UNA vez y se reutilizan para HO
%   y para todos los p de HRS (son subconjuntos de los mismos niveles).
[~, N] = size(V);
ps = cfg.hurst.particiones_hrs;
mmin = cfg.hurst.bloque_min_abs;
todos = unique(floor(N ./ 2.^(0:floor(log2(N)))), 'stable');
todos = todos(todos >= mmin);
RS = rs_diadico(V, todos, cfg.hurst.promedio_rs);

ms_ho = todos(todos >= max(cfg.hurst.ho_bloque_min, mmin));
[Hho, nho] = pendiente_hurst(ms_ho, RS(:, ismember(todos, ms_ho)), cfg.hurst.min_escalas);

Hhrs = zeros(size(V,1), numel(ps)); nhrs = Hhrs;
for i = 1:numel(ps)
    ms = niveles_rangors(N, ps(i), mmin);
    [Hhrs(:,i), nhrs(:,i)] = pendiente_hurst(ms, RS(:, ismember(todos, ms)), cfg.hurst.min_escalas);
end
[Hhv, nhv] = hurst_hv(V, cfg);

f    = [Hho.', Hhrs(:).', Hhv.'];
nesc = [nho.', nhrs(:).', nhv.'];
end
