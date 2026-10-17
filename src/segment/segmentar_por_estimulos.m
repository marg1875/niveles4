function A = segmentar_por_estimulos(N, stims, cfg)
%SEGMENTAR_POR_ESTIMULOS  Bloques y anclas temporales a partir de stims.
%   A = segmentar_por_estimulos(N, stims, cfg) devuelve una tabla con una
%   fila por ancla: bloque, tarea (false = reposo, true = tarea), ancla
%   (1..3) y centro (muestra, base 1). El NIVEL (10/40) no se toma de stims
%   (en algunas sesiones el programa no lo actualizó) sino del nombre.
%
%   stims NO trae marcas temporales: es la secuencia de etiquetas mostradas
%   ([0 0 0 0 L L L L] x 10). Se asume espaciado uniforme (respaldado por
%   qa_alineacion_estimulos.m): la entrada k ocupa las muestras
%   [round((k-1)N/K)+1, round(kN/K)]. Las corridas de etiquetas iguales
%   forman bloques (~7.8 s). En cada bloque se descartan márgenes
%   (transición/reacción) y se fijan cfg.segmentacion.anclas_por_bloque
%   anclas: centros de ranuras contiguas de 512 muestras centradas en la
%   zona útil. Toda duración de ventana se centra en estas mismas anclas,
%   así que todas las duraciones analizan los mismos instantes y el mismo
%   número de instancias.
fs = cfg.fs;
K = numel(stims);
bordes = round((0:K) * N / K);              % la entrada k: bordes(k)+1 .. bordes(k+1)
ini = [1, find(diff(stims) ~= 0) + 1];      % primera entrada de cada bloque
fin = [ini(2:end) - 1, K];
m0 = round(cfg.segmentacion.margen_inicio_s * fs);
m1 = round(cfg.segmentacion.margen_fin_s * fs);
R = cfg.segmentacion.muestras_ranura;

bloque = []; tarea = []; ancla = []; centro = [];
for b = 1:numel(ini)
    u0 = bordes(ini(b)) + 1 + m0;
    u1 = bordes(fin(b) + 1) - m1;
    largo = u1 - u0 + 1;
    na = min(cfg.segmentacion.anclas_por_bloque, floor(largo / R));
    if na < 1, continue, end
    desplaz = u0 + floor((largo - na*R) / 2);
    for a = 1:na
        bloque(end+1,1) = b; %#ok<AGROW>
        tarea(end+1,1)  = stims(ini(b)) > 0; %#ok<AGROW>
        ancla(end+1,1)  = a; %#ok<AGROW>
        centro(end+1,1) = desplaz + (a-1)*R + R/2; %#ok<AGROW>
    end
end
A = table(bloque, logical(tarea), ancla, centro, 'VariableNames', {'bloque','tarea','ancla','centro'});
end
