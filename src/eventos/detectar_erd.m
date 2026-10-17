function E = detectar_erd(X, cfg)
%DETECTAR_ERD  Detección de intentos motores por desincronización (ERD).
%   E = detectar_erd(X, cfg), con X (16 x N) la señal limpia de un archivo.
%   Pensado para la sesión imaginada (KMI), que fue ASÍNCRONA: el paciente
%   intentó el ejercicio cuantas veces quiso, así que stims no marca cuándo.
%   Port de Niveles3/src/events/detector.py:
%     1) potencia mu/beta (cfg.eventos.banda, Butterworth de orden 4,
%        filtfilt) promediada en cfg.eventos.canales_erd y suavizada con una
%        media móvil de 1 s;
%     2) umbral autocalibrado = p10 + 0.2 (p90 - p10);
%     3) evento = tramo por debajo del umbral de al menos
%        cfg.eventos.min_duracion_s; se ignoran los primeros
%        cfg.eventos.reposo_inicial_s (reposo inicial) y se fusionan/
%        descartan eventos cuyos mínimos estén a menos de
%        cfg.eventos.separacion_min_s (se conserva el más profundo).
%   E: tabla con inicio, fin y minimo (muestras, base 1) y profundidad
%   (energía relativa al umbral).
fs = cfg.fs;
idx = find(ismember(cfg.canales.nombres, cfg.eventos.canales_erd));
[b, a] = butter(4, cfg.eventos.banda / (fs/2), 'bandpass');
Y = filtfilt(b, a, double(X(idx, :)).');
P = mean(Y.^2, 2).';
env = movmean(P, fs);
q = prctile_simple(env, [10 90]);
umbral = q(1) + cfg.eventos.factor_umbral * (q(2) - q(1));
bajo = env < umbral;
bajo(1:round(cfg.eventos.reposo_inicial_s * fs)) = false;
d = diff([false, bajo, false]);
ini = find(d == 1);
fin = find(d == -1) - 1;
largo = fin - ini + 1;
ok = largo >= round(cfg.eventos.min_duracion_s * fs);
ini = ini(ok); fin = fin(ok);
minimo = zeros(size(ini)); prof = zeros(size(ini));
for k = 1:numel(ini)
    [v, j] = min(env(ini(k):fin(k)));
    minimo(k) = ini(k) + j - 1;
    prof(k) = v / umbral;
end
% Separación mínima entre eventos: se conserva el más profundo
sep = round(cfg.eventos.separacion_min_s * fs);
[~, orden] = sort(prof);
quedan = false(size(ini));
for k = orden
    if ~any(quedan & abs(minimo - minimo(k)) < sep)
        quedan(k) = true;
    end
end
E = table(ini(quedan).', fin(quedan).', minimo(quedan).', prof(quedan).', ...
    'VariableNames', {'inicio','fin','minimo','profundidad'});
E = sortrows(E, 'minimo');
end

function q = prctile_simple(x, p)
% Percentiles sin Statistics Toolbox (interpolación lineal).
x = sort(x(:));
n = numel(x);
pos = 1 + (n - 1) * p / 100;
q = interp1(1:n, x, pos);
end
