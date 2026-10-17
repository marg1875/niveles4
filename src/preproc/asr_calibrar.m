function E = asr_calibrar(X, cfg)
%ASR_CALIBRAR  Estado de calibración de Artifact Subspace Reconstruction.
%   E = asr_calibrar(X, cfg), X (canales x muestras) ya filtrada.
%   Sigue a Mullen et al. (2015, IEEE TBME 62:2553) / clean_rawdata
%   (asr_calibrate) con dos simplificaciones documentadas:
%     - Sin prefiltro espectral (la señal ya viene en 0.5-45 Hz).
%     - Las distribuciones de RMS se resumen con mediana y 1.4826*MAD en
%       lugar del ajuste de gaussiana truncada.
%   1) Datos de calibración = ventanas "limpias" del propio registro: RMS
%      por canal en ventanas de cfg.asr.ventana_s; una ventana es limpia si
%      en >= cfg.asr.frac_canales de los canales su z robusto cae dentro de
%      cfg.asr.z_calib.
%   2) C = mediana (elemento a elemento) de las covarianzas de esas
%      ventanas; M = sqrtm(C); V = autovectores de M.
%   3) Umbral por componente: T = diag(mu + k*sigma) * V', con mu/sigma del
%      RMS de las componentes V'*x en las ventanas limpias.
%   E: struct con M, T, k, frac_limpia (fracción de ventanas de calibración).
N = max(2, round(cfg.asr.ventana_s * cfg.fs));
[nc, ns] = size(X);
X = X - median(X, 2);
ini = 1:N:ns-N+1;
nv = numel(ini);
rms = zeros(nc, nv);
for v = 1:nv
    rms(:, v) = sqrt(mean(X(:, ini(v):ini(v)+N-1).^2, 2));
end
med = median(rms, 2);
z = (rms - med) ./ max(1.4826 * median(abs(rms - med), 2), eps);
ok = mean(z >= cfg.asr.z_calib(1) & z <= cfg.asr.z_calib(2), 1) >= cfg.asr.frac_canales;
if nnz(ok) < max(nc, 10)                 % registro muy sucio: se usa todo
    ok = true(1, nv);
end
vs = find(ok);
covs = zeros(nc, nc, numel(vs));
for i = 1:numel(vs)
    B = X(:, ini(vs(i)):ini(vs(i))+N-1);
    covs(:, :, i) = (B * B.') / N;
end
C = median(covs, 3);
C = (C + C.') / 2;
M = real(sqrtm(C));
[V, D] = eig((M + M.') / 2);
[~, orden] = sort(diag(D));
V = V(:, orden);
% RMS de cada componente en las ventanas limpias
rc = zeros(nc, numel(vs));
for i = 1:numel(vs)
    B = V.' * X(:, ini(vs(i)):ini(vs(i))+N-1);
    rc(:, i) = sqrt(mean(B.^2, 2));
end
mu = median(rc, 2);
sg = 1.4826 * median(abs(rc - mu), 2);
E.M = M;
E.T = diag(mu + cfg.asr.k * sg) * V.';
E.k = cfg.asr.k;
E.frac_limpia = mean(ok);
end
