function etapa2_caracteristicas(cfg)
%ETAPA2_CARACTERISTICAS  Segmentación multi-duración + HO/HRS(p)/HV.
%   etapa2_caracteristicas(cfg) lee las señales limpias de la etapa 1 y, para
%   cada duración L de cfg.segmentacion.duraciones, extrae una ventana de L
%   muestras centrada en cada ancla (3 por bloque de reposo/tarea) y calcula
%   por canal: HO, HRS(p) para cada p, y HV.
%   -> resultados/caracteristicas/caracteristicas_d<L>.csv
%      (una fila por ventana: metadatos + 144 características)
%   -> resultados/caracteristicas/escalas_d<L>.csv (n.º de escalas usadas)

dirS = fullfile(cfg.rutas.resultados, 'senales');
dirC = fullfile(cfg.rutas.resultados, 'caracteristicas');
if ~isfolder(dirC), mkdir(dirC); end
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');
if cfg.modo_prueba
    ses = ses(strcmp(ses.paciente, 'Px006') & ses.mes == 1, :);
end
duraciones = cfg.segmentacion.duraciones;
pend = duraciones(~arrayfun(@(L) isfile(fullfile(dirC, sprintf('caracteristicas_d%d.csv', L))), duraciones));
if isempty(pend)
    registrar(cfg, 'Etapa 2: características ya calculadas para todas las duraciones.');
    return
end
registrar(cfg, 'Etapa 2: %d sesiones x %d duraciones.', height(ses), numel(pend));
nombres = nombres_caracteristicas(cfg);
nw = usar_paralelo(cfg);
partes = cell(height(ses), 1);
sesion = ses.sesion; paciente = ses.paciente; mes = ses.mes; condicion = ses.condicion; completa = ses.completa;
parfor (s = 1:height(ses), nw)
    S = load(fullfile(dirS, [sesion{s} '.mat']));
    salida = cell(1, numel(pend));
    for j = 1:numel(S.X)
        [~, arch] = fileparts(S.archivos{j});
        A0 = segmentar_por_estimulos(size(S.X{j}, 2), S.stims{j}, cfg);
        malas = zeros(0, 2);
        if isfield(S, 'malas'), malas = S.malas{j}; end
        for d = 1:numel(pend)
            L = pend(d);
            % rechazo de artefactos: se descartan las ventanas que tocan tramos malos
            A = A0(ventana_valida(malas, A0.centro - floor(L/2), L), :);
            clase = A.tarea * S.niveles(j);              % 0 = reposo, 10/40 = tarea
            Vs = extraer_ventanas(S.X{j}, A.centro, L);
            nv = height(A);
            F = zeros(nv, numel(nombres)); NE = F;
            for e = 1:nv
                [F(e,:), NE(e,:)] = extraer_hurst_ventana(double(Vs(:,:,e)), cfg);
            end
            meta = table(repmat(sesion(s), nv, 1), repmat(paciente(s), nv, 1), ...
                repmat(mes(s), nv, 1), repmat(condicion(s), nv, 1), repmat(completa(s), nv, 1), ...
                repmat({arch}, nv, 1), repmat(S.niveles(j), nv, 1), clase, A.bloque, A.ancla, ...
                repmat(L, nv, 1), 'VariableNames', {'sesion','paciente','mes','condicion', ...
                'sesion_completa','archivo','nivel_archivo','clase','bloque','ancla','duracion'});
            salida{d} = [salida{d}; [meta, array2table(F, 'VariableNames', nombres), ...
                array2table(NE, 'VariableNames', strcat('n_', nombres))]];
        end
    end
    partes{s} = salida;
    fprintf('  %s listo\n', sesion{s});
end
for d = 1:numel(pend)
    T = vertcat(cellfun(@(c) c{d}, partes, 'UniformOutput', false));
    T = vertcat(T{:});
    esc = T(:, [1:11, find(startsWith(T.Properties.VariableNames, 'n_'))]);
    T(:, startsWith(T.Properties.VariableNames, 'n_')) = [];
    writetable(esc, fullfile(dirC, sprintf('escalas_d%d.csv', pend(d))), 'Encoding', 'UTF-8');
    writetable(T, fullfile(dirC, sprintf('caracteristicas_d%d.csv', pend(d))), 'Encoding', 'UTF-8');
    fn = T{:, nombres};
    registrar(cfg, '  d=%d muestras (%.2f s): %d ventanas; NaN = %.1f %%; HO medio = %.3f, HV medio = %.3f', ...
        pend(d), pend(d)/cfg.fs, height(T), 100*mean(isnan(fn(:))), ...
        mean(fn(:, startsWith(nombres,'HO_')), 'all', 'omitnan'), ...
        mean(fn(:, startsWith(nombres,'HV_')), 'all', 'omitnan'));
end
end
