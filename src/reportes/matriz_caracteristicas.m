function X = matriz_caracteristicas(T)
%MATRIZ_CARACTERISTICAS  Columna 'feat' de una tabla de ventanas como matriz.
if iscell(T.feat), X = vertcat(T.feat{:}); else, X = T.feat; end
end
