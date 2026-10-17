function F = disenar_filtros(cfg)
%DISENAR_FILTROS  Diseña una sola vez los filtros del preprocesamiento.
%   F = disenar_filtros(cfg) devuelve una estructura con:
%     F.notch_b, F.notch_a : notch IIR de 2.º orden a cfg.filtro.notch_hz
%                            (misma fórmula que scipy.signal.iirnotch, usada
%                            en Niveles3/src/preprocessing/filters.py)
%     F.tipo = 'butter' -> F.sos, F.g  (Butterworth pasabanda en SOS)
%     F.tipo = 'fir'    -> F.b          (FIR fir1 pasabanda)
fs = cfg.fs;

% --- Notch (implementado a mano para no depender de DSP System Toolbox)
w0 = 2*pi*cfg.filtro.notch_hz/fs;
bw = w0 / cfg.filtro.notch_q;
beta = tan(bw/2);                 % ganancia -3 dB en los bordes de la banda
ganancia = 1/(1+beta);
F.notch_b = ganancia * [1, -2*cos(w0), 1];
F.notch_a = [1, -2*ganancia*cos(w0), 2*ganancia-1];

% --- Pasabanda
F.tipo = cfg.filtro.tipo;
Wn = cfg.filtro.banda/(fs/2);
switch cfg.filtro.tipo
    case 'butter'
        % SOS: el pasabanda de orden 4 con corte en 0.5 Hz a 250 Hz está mal
        % condicionado en forma de función de transferencia.
        [z,p,k] = butter(cfg.filtro.orden_butter, Wn, 'bandpass');
        [F.sos, F.g] = zp2sos(z,p,k);
    case 'fir'
        F.b = fir1(cfg.filtro.orden_fir, Wn, 'bandpass');
    otherwise
        error('disenar_filtros:tipo', 'Tipo de filtro desconocido: %s', cfg.filtro.tipo);
end
end
