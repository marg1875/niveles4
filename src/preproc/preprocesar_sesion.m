function S = preprocesar_sesion(rutas, cfg, F)
%PREPROCESAR_SESION  Preprocesa los archivos de una sesión con una única W.
%   S = preprocesar_sesion(rutas, cfg, F), con rutas = celda de .mat de la
%   misma sesión (paciente x mes x condición: archivo de 10 % y de 40 %).
%
%   1) Por archivo: DC -> notch -> pasabanda -> puentes (cfg.puentes) ->
%      ASR (cfg.asr) -> CAR (cfg.preproc.car) ->
%      z-score por canal (Ec. 1).
%   4) Se marcan los tramos con artefactos residuales (marcar_artefactos):
%      S.malas{j} = [inicio fin]; las etapas siguientes descartan las
%      ventanas que los tocan (ventana_valida).
%   2) Por sesión (archivos concatenados): blanqueo (Ec. 2) + ICA (Ecs. 4-6).
%      Una sola W por sesión hace que las tres clases (reposo/10/40) pasen
%      por el mismo desmezclador: con una W por archivo, el clasificador
%      podría aprender "qué archivo" en lugar de "qué nivel".
%   3) Se anulan las fuentes artefactuales y se reproyecta cada archivo al
%      espacio de electrodos: X_limpio = X - A(:,idx)*Y(idx,:), A = pinv(W*Q).
%
%   Devuelve S.X{j} (16 x N_j, single), S.stims{j}, S.info.
nf = numel(rutas);
Xs = cell(1, nf);
Xf = cell(1, nf);
S.stims = cell(1, nf);
for j = 1:nf
    [X, st] = cargar_registro(rutas{j});
    Xf{j} = filtrar_eeg(X, F);
    S.stims{j} = st;
end
% Electrodos puenteados (cfg.puentes): se detectan por archivo; con
% 'eliminar' se usa la unión de la sesión para que el ICA tenga las mismas
% dimensiones en todos los archivos.
metodo = 'ninguno';
if isfield(cfg, 'puentes'), metodo = cfg.puentes.metodo; end
S.info.puenteados = cell(1, nf);
if ~strcmp(metodo, 'ninguno')
    c0 = cfg; c0.puentes.metodo = 'ninguno';       % solo detectar
    for j = 1:nf
        [~, S.info.puenteados{j}] = corregir_puentes(Xf{j}, c0);
    end
    for j = 1:nf
        p = S.info.puenteados{j};
        if strcmp(metodo, 'eliminar'), p = unique([S.info.puenteados{:}]); end
        Xf{j} = corregir_puentes(Xf{j}, cfg, p);
    end
end
% ASR (cfg.asr): repara ráfagas de alta varianza por archivo, antes del CAR
S.info.frac_asr = zeros(1, nf);
if isfield(cfg, 'asr') && cfg.asr.activar
    for j = 1:nf
        ok = all(isfinite(Xf{j}), 2);              % canales eliminados por puenteo quedan fuera
        E = asr_calibrar(Xf{j}(ok, :), cfg);
        [Xf{j}(ok, :), S.info.frac_asr(j)] = asr_procesar(Xf{j}(ok, :), E, cfg);
    end
end
for j = 1:nf
    if isfield(cfg, 'preproc') && cfg.preproc.car
        Xf{j} = referencia_promedio(Xf{j});       % CAR antes de normalizar y del ICA
    end
    Xs{j} = normalizar_zscore(Xf{j});
end
validos = all(cell2mat(cellfun(@(x) all(isfinite(x), 2), Xs, 'UniformOutput', false)), 2);
S.info.ica = cfg.ica.activar;
S.info.removidos = [];
S.info.curtosis = [];
S.info.iteraciones = 0;
S.info.convergio = true;
S.info.rango = nnz(validos);

if cfg.ica.activar && nnz(validos) < 3
    % casi todos los canales eliminados por puenteo: sin ICA (la sesión
    % queda excluida por la regla de calidad en los análisis)
    S.info.ica = false;
elseif cfg.ica.activar
    Xc = cell2mat(cellfun(@(x) x(validos, :), Xs, 'UniformOutput', false));
    % W se estima con una de cada cfg.ica.submuestreo muestras (la señal ya
    % está limitada a 45 Hz y el modelo ICA es instantáneo: el orden temporal
    % no interviene) y luego se aplica a TODAS las muestras.
    Xe = Xc(:, 1:cfg.ica.submuestreo:end);
    [~, Q, binfo] = blanquear(Xe, cfg);
    Z = Q * (Xe - mean(Xe, 2));
    [W, ~, iinfo] = ica_gradiente_natural(Z, cfg);
    Y = (W * Q) * (Xc - mean(Xe, 2));
    [idx, curt] = limpiar_artefactos_ica(Y, cfg);
    A = pinv(W * Q);
    Wt = W * Q;                                    % fuentes desde datos z-score
    for j = 1:nf
        if ~isempty(idx)                           % solo canales válidos (no eliminados)
            Yj = Wt(idx,:) * (Xs{j}(validos, :) - mean(Xe, 2));
            Xs{j}(validos, :) = Xs{j}(validos, :) - A(:, idx) * Yj;
        end
    end
    S.info.removidos = idx;
    S.info.curtosis = curt.';
    S.info.iteraciones = iinfo.iteraciones;
    S.info.convergio = iinfo.convergio;
    S.info.rango = binfo.rango;
end
S.X = cellfun(@single, Xs, 'UniformOutput', false);
% Tramos con artefactos residuales (tras CAR + ICA), por archivo
S.malas = cell(1, nf);
S.info.frac_mala = zeros(1, nf);
for j = 1:nf
    if isfield(cfg, 'artefactos')
        S.malas{j} = marcar_artefactos(Xs{j}, cfg);
    else
        S.malas{j} = zeros(0, 2);
    end
    S.info.frac_mala(j) = sum(diff(S.malas{j}, 1, 2) + 1) / size(Xs{j}, 2);
end
end
