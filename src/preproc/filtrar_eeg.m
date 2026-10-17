function Y = filtrar_eeg(X, F)
%FILTRAR_EEG  Remoción de DC, notch y pasabanda de fase cero (filtfilt).
%   Y = filtrar_eeg(X, F) con X de tamaño canales x muestras y F de
%   disenar_filtros. Se filtra el registro CONTINUO completo (antes de
%   segmentar) para evitar efectos de borde en las ventanas cortas.
Xt = X.' - mean(X, 2).';                       % muestras x canales, sin DC
Xt = filtfilt(F.notch_b, F.notch_a, Xt);
switch F.tipo
    case 'butter', Xt = filtfilt(F.sos, F.g, Xt);
    case 'fir',    Xt = filtfilt(F.b, 1, Xt);
end
Y = Xt.';
end
