function E = eventos_fis(N, stims)
%EVENTOS_FIS  Bloques de estímulo de la sesión física (síncrona).
%   E = eventos_fis(N, stims) devuelve una tabla con inicio, fin y minimo
%   (aquí: centro del bloque) de cada bloque de TAREA (stims > 0), con la
%   misma correspondencia uniforme entrada->muestras que
%   segmentar_por_estimulos: la entrada k ocupa [round((k-1)N/K)+1, round(kN/K)].
K = numel(stims);
bordes = round((0:K) * N / K);
tarea = stims(:).' > 0;
d = diff([false, tarea, false]);
ki = find(d == 1);
kf = find(d == -1) - 1;
inicio = bordes(ki) + 1;
fin = bordes(kf + 1);
minimo = round((inicio + fin) / 2);
E = table(inicio(:), fin(:), minimo(:), nan(numel(inicio), 1), ...
    'VariableNames', {'inicio','fin','minimo','profundidad'});
end
