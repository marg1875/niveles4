function [X, puenteados] = corregir_puentes(X, cfg, puenteados)
%CORREGIR_PUENTES  Detecta y corrige electrodos puenteados (gel entre electrodos).
%   [X, puenteados] = corregir_puentes(X, cfg) con X (canales x muestras),
%   ya filtrada y ANTES de la referencia promedio. Dos canales están
%   puenteados si, tras aplicar CAR, su correlación supera
%   cfg.puentes.umbral_r (en EEG genuino dos electrodos distintos no
%   registran señales casi idénticas una vez quitado el modo común).
%   cfg.puentes.metodo:
%     'ninguno'    -> no se toca la señal
%     'eliminar'   -> los canales puenteados quedan en NaN (se excluyen del
%                     ICA y sus características quedan como faltantes '?')
%     'interpolar' -> cada canal puenteado se reemplaza por el promedio,
%                     ponderado por la inversa de la distancia, de sus
%                     cfg.puentes.vecinos vecinos NO puenteados más cercanos
%   [X, p] = corregir_puentes(X, cfg, p) aplica una lista de canales ya
%   detectada (útil para usar la misma lista en todos los archivos de una sesión).
if nargin < 3
    Y = X - mean(X, 1);
    R = corrcoef(Y.');
    R(1:size(R, 1) + 1:end) = 0;
    puenteados = find(any(abs(R) > cfg.puentes.umbral_r, 2)).';
end
if isempty(puenteados) || strcmp(cfg.puentes.metodo, 'ninguno'), return, end
switch cfg.puentes.metodo
    case 'eliminar'
        X(puenteados, :) = NaN;
    case 'interpolar'
        P = coordenadas_canales(cfg.canales.nombres);
        buenos = setdiff(1:size(X, 1), puenteados);
        X0 = X;
        for c = puenteados
            d = sqrt(sum((P(buenos, :) - P(c, :)).^2, 2));
            [d, o] = sort(d);
            k = min(cfg.puentes.vecinos, numel(o));
            w = 1 ./ max(d(1:k), eps);
            X(c, :) = (w.' * X0(buenos(o(1:k)), :)) / sum(w);
        end
end
end
