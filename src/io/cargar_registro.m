function [X, stims] = cargar_registro(ruta)
%CARGAR_REGISTRO  Carga un registro (formato A o legado).
%   [X, stims] = cargar_registro(ruta) devuelve X (16 x N, double) y stims.
%   Formato A: X = EEG_all y stims (1 x 80). La variable B se ignora: es una
%   copia de las últimas 2500 muestras de EEG_all (verificado).
%   Formato legado (Px009 mes 1 y px010_10_1; solo B, ~89 s): X = B y
%   stims = [] porque su secuencia de estímulos (0 L 0 L ..., 320 entradas
%   en KMI) no se puede alinear con la señal. Estos registros entran a los
%   análisis que usan la etiqueta del archivo (10 vs 40) y no a los que
%   necesitan los bloques de reposo/tarea (segmentar_por_estimulos los omite).
S = load(ruta);
if isfield(S, 'EEG_all')
    X = double(S.EEG_all);
    stims = double(S.stims(:)');
elseif isfield(S, 'B')
    X = double(S.B);
    stims = [];
else
    error('cargar_registro:formato', '%s no contiene EEG_all ni B.', ruta);
end
X(~isfinite(X)) = 0;
end
