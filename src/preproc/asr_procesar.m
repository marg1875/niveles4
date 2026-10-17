function [Y, frac] = asr_procesar(X, E, cfg)
%ASR_PROCESAR  Repara ráfagas de alta varianza con Artifact Subspace Reconstruction.
%   [Y, frac] = asr_procesar(X, E, cfg), X (canales x muestras), E de
%   asr_calibrar. Ventanas de cfg.asr.ventana_s centradas cada media ventana:
%     [V, D] = eig(cov(ventana)); una componente se conserva si su varianza
%     D(i) < sum((T*V(:,i)).^2) (umbral de calibración proyectado), y nunca se
%     reparan más de cfg.asr.max_dims*canales componentes (las de mayor
%     varianza). R = M * pinv(diag(keep) * V' * M) * V' reconstruye el
%     subespacio rechazado a partir del resto (Mullen et al., 2015).
%   Entre dos centros, la salida se mezcla con coseno elevado de R_k a R_k+1
%   (fuera de las ráfagas R = I y la señal queda intacta). Versión no causal
%   (registros completos, sin retardo).
%   frac: fracción de muestras afectadas por alguna reconstrucción.
N = max(2, round(cfg.asr.ventana_s * cfg.fs));
paso = max(1, round(N / 2));
[nc, ns] = size(X);
centros = 1:paso:ns;
if centros(end) < ns, centros(end+1) = ns; end
maxd = round(cfg.asr.max_dims * nc);
nk = numel(centros);
R = cell(1, nk);
triv = true(1, nk);
for k = 1:nk
    b = min(ns, max(1, centros(k) - floor(N/2)) + N - 1);
    a = max(1, b - N + 1);
    B = X(:, a:b);
    Cw = (B * B.') / size(B, 2);
    [V, D] = eig((Cw + Cw.') / 2);
    [d, o] = sort(diag(D));
    V = V(:, o);
    keep = d.' < sum((E.T * V).^2, 1) | (1:nc) <= nc - maxd;
    if ~all(keep)
        R{k} = real(E.M * pinv(keep.' .* (V.' * E.M)) * V.');
        triv(k) = false;
    end
end
Y = X;
afectada = false(1, ns);
I = eye(nc);
for k = 1:nk-1
    if triv(k) && triv(k+1), continue, end
    t = centros(k):centros(k+1);
    w = (1 - cos(pi * (t - centros(k)) / (centros(k+1) - centros(k)))) / 2;
    Ra = I; if ~triv(k), Ra = R{k}; end
    Rb = I; if ~triv(k+1), Rb = R{k+1}; end
    Y(:, t) = (1 - w) .* (Ra * X(:, t)) + w .* (Rb * X(:, t));
    afectada(t) = true;
end
frac = mean(afectada);
end
