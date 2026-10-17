function [X, stims] = cargar_registro(ruta)
%CARGAR_REGISTRO  Carga un registro de formato A.
%   [X, stims] = cargar_registro(ruta) devuelve X = EEG_all (16 x N, double)
%   y stims (1 x 80, double). La variable B se ignora: es una copia de las
%   últimas 2500 muestras de EEG_all (verificado en la exploración).
S = load(ruta, 'EEG_all', 'stims');
if ~isfield(S, 'EEG_all')
    error('cargar_registro:formato', '%s no contiene EEG_all (formato legado).', ruta);
end
X = double(S.EEG_all);
X(~isfinite(X)) = 0;
stims = double(S.stims(:)');
end
