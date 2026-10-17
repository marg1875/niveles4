function [W, Y, info] = ica_gradiente_natural(Z, cfg)
%ICA_GRADIENTE_NATURAL  ICA por gradiente natural con no linealidad flexible.
%   [W, Y, info] = ica_gradiente_natural(Z, cfg) estima la matriz de
%   separación W (k x k) para datos blanqueados Z (k x n): Y = W*Z.
%
%   Algoritmo (Ecs. 5-6 del artículo; Cichocki & Amari 2002):
%     W <- W + eta * [Lambda - E{phi(y) y'}] * W
%   con la función de activación flexible según la curtosis kappa_i:
%     kappa_i <= 0 (sub-gaussiana):  phi = alfa*lam*sgn(y)*|lam*y|^(alfa-1)
%     kappa_i  > 0 (super-gaussiana): phi = (1+beta)*y / (y^2 + beta/lam^2)
%   con lam_i = 1/sigma(y_i).
%
%   Detalles de implementación (documentados en el README):
%   - Se usa la versión "no holonómica" Lambda = diag(E{phi_i(y_i) y_i}):
%     solo los términos cruzados (dependencia entre componentes) impulsan el
%     aprendizaje; la escala de cada fuente es indeterminada en ICA, así que
%     tras cada paso se renormaliza cada fila de W a varianza unitaria.
%   - Cada fila del gradiente se escala por 1/Lambda_ii para que eta sea
%     adimensional (independiente de alfa/beta).
%   - W0 = I (determinista y reproducible, a diferencia de FastICA).
%   - Paso adaptativo: crece 5 % mientras max|gradiente| disminuya (tope
%     cfg.ica.eta_max) y se reduce 30 % si aumenta (piso cfg.ica.eta).
%   - Paro cuando max|gradiente relativo| < cfg.ica.tol.
%   - Si diverge (NaN/Inf), se reinicia con eta/2 (cfg.ica.max_reintentos).

k = size(Z, 1);
n = size(Z, 2);
eta = cfg.ica.eta;
for intento = 0:cfg.ica.max_reintentos
    W = eye(k);
    convergio = false;
    diverge = false;
    historial = nan(cfg.ica.max_iter, 1);
    paso = eta;
    g_prev = Inf;
    for it = 1:cfg.ica.max_iter
        Y = W * Z;
        sd = std(Y, 0, 2);
        lam = 1 ./ max(sd, eps);
        kap = mean(Y.^4, 2) ./ max(mean(Y.^2, 2).^2, eps) - 3;
        phi = zeros(size(Y));
        sub = kap <= 0;
        if any(sub)
            U = lam(sub) .* Y(sub,:);
            phi(sub,:) = cfg.ica.alfa .* lam(sub) .* sign(U) .* abs(U).^(cfg.ica.alfa-1);
        end
        if any(~sub)
            Ys = Y(~sub,:);
            phi(~sub,:) = (1+cfg.ica.beta) .* Ys ./ (Ys.^2 + cfg.ica.beta ./ lam(~sub).^2);
        end
        R = (phi * Y.') / n;                       % E{phi(y) y'}
        Ld = diag(R);
        G = diag(Ld) - R;                          % regla no holonómica
        G = G ./ max(abs(Ld), eps);                % escala por fila
        dW = paso * G * W;
        W = W + dW;
        W = W ./ std(W*Z, 0, 2);                   % varianza unitaria por fuente
        if any(~isfinite(W(:)))
            diverge = true; break
        end
        % Criterio de paro sobre el GRADIENTE (independiente del paso): en el
        % punto fijo E{phi_i(y_i) y_j} = 0 para i ~= j.
        g = max(abs(G(:)));
        historial(it) = g;
        % Paso adaptativo ("bold driver") acotado a [eta, eta_max]: con EEG de
        % colas pesadas la convergencia con paso fijo es lineal y muy lenta.
        if g < g_prev
            paso = min(paso * 1.05, cfg.ica.eta_max);
        else
            paso = max(paso * 0.7, eta);
        end
        g_prev = g;
        if g < cfg.ica.tol
            convergio = true; break
        end
    end
    if ~diverge, break, end
    eta = eta/2;
end
if diverge
    warning('ica_gradiente_natural:diverge', 'ICA divergió; se devuelve W = I.');
    W = eye(k);
end
Y = W * Z;
info.iteraciones = it;
info.historial = historial(1:it);
info.convergio = convergio;
info.eta_final = eta;
info.curtosis = mean(Y.^4, 2) ./ mean(Y.^2, 2).^2 - 3;
end
