function W = ventanas_por_definicion(E, N, def, condicion, cfg)
%VENTANAS_POR_DEFINICION  Centros de ventanas de tarea y reposo de un archivo.
%   W = ventanas_por_definicion(E, N, def, condicion, cfg)
%     E   : eventos del archivo (eventos_fis o detectar_erd)
%     def : 'A' (reposo = primeros segundos del registro) o
%           'B' (reposo = [fin + 0.5 s, fin + 0.8 s] tras cada evento)
%   Todas las ventanas de un mismo dataset tienen la MISMA longitud L
%   (cfg.eventos.L_A o cfg.eventos.L_B), porque el exponente de Hurst
%   depende de la longitud de la ventana.
%   Tarea: FIS -> cfg.eventos.anclas_fis ventanas por bloque, centradas en
%   ranuras equiespaciadas dentro del bloque (tras un margen de 0.5 s);
%   KMI -> una ventana centrada en el mínimo de energía de cada ERD.
%   Devuelve una tabla: tipo ('tarea'|'reposo'), evento, inicio (muestra
%   de inicio de la ventana, base 1) y L.
fs = cfg.fs;
switch def
    case 'A', L = cfg.eventos.L_A;
    case 'B', L = cfg.eventos.L_B;
end
tipo = {}; evento = []; inicio = [];
% --- Tarea
for k = 1:height(E)
    if strcmp(condicion, 'FIS')
        u0 = E.inicio(k) + round(0.5 * fs);
        u1 = E.fin(k);
        na = cfg.eventos.anclas_fis;
        c = round(u0 + ((1:na) - 0.5) * (u1 - u0) / na);
    else
        c = E.minimo(k);
    end
    s = c - floor(L/2);
    s = s(s >= 1 & s + L - 1 <= N);
    tipo = [tipo, repmat({'tarea'}, 1, numel(s))]; %#ok<AGROW>
    evento = [evento, repmat(k, 1, numel(s))]; %#ok<AGROW>
    inicio = [inicio, s]; %#ok<AGROW>
end
% --- Reposo
switch def
    case 'A'
        s = ventanas_reposo_inicial(L, cfg);
        s = s(s + L - 1 <= N);
        tipo = [tipo, repmat({'reposo'}, 1, numel(s))];
        evento = [evento, zeros(1, numel(s))];
        inicio = [inicio, s];
    case 'B'
        for k = 1:height(E)
            s = E.fin(k) + round(cfg.eventos.post_ini_s * fs);
            if s + L - 1 <= N && (k == height(E) || s + L - 1 < E.inicio(k+1))
                tipo{end+1} = 'reposo'; evento(end+1) = k; inicio(end+1) = s; %#ok<AGROW>
            end
        end
end
W = table(tipo(:), evento(:), inicio(:), repmat(L, numel(inicio), 1), ...
    'VariableNames', {'tipo','evento','inicio','L'});
end
