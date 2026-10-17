function sesiones = etapa1_preprocesar(cfg)
%ETAPA1_PREPROCESAR  Inventario, carga y preprocesamiento por sesión.
%   sesiones = etapa1_preprocesar(cfg)
%   1) Inventario de Datos/ -> resultados/registros.csv (con exclusiones y
%      avisos).
%   2) Agrupa los archivos válidos en SESIONES (paciente x mes x condición).
%      Una sesión es "completa" si tiene el archivo de 10 % y el de 40 %.
%   3) Preprocesa cada sesión con una única W (preprocesar_sesion) y guarda
%      las señales limpias en resultados/senales/<sesion>.mat, de modo que
%      el ICA se calcula una sola vez para todos los barridos posteriores.
%   -> resultados/sesiones.csv

dirS = fullfile(cfg.rutas.resultados, 'senales');
if ~isfolder(dirS), mkdir(dirS); end

T = listar_registros(cfg);
writetable(T(:, setdiff(T.Properties.VariableNames, {'ruta'}, 'stable')), ...
    fullfile(cfg.rutas.resultados, 'registros.csv'), 'Encoding', 'UTF-8');
registrar(cfg, 'Inventario: %d archivos, %d válidos, %d excluidos.', height(T), ...
    nnz(strcmp(T.estado,'ok')), nnz(~strcmp(T.estado,'ok')));

V = T(strcmp(T.estado, 'ok'), :);
if cfg.modo_prueba
    V = V(strcmp(V.paciente, 'Px006') & V.mes == 1, :);
end
[g, pac, mes, cond] = findgroups(V.paciente, V.mes, V.condicion);
ns = max(g);
id = cell(ns,1); completa = false(ns,1); archivos = cell(ns,1); niveles = cell(ns,1);
for s = 1:ns
    id{s} = sprintf('%s_M%d_%s', pac{s}, mes(s), cond{s});
    filas = V(g == s, :);
    filas = sortrows(filas, 'nivel');
    archivos{s} = filas.ruta.';
    niveles{s} = filas.nivel.';
    completa(s) = all(ismember([10 40], filas.nivel));
end
sesiones = table(id, pac, mes, cond, completa, archivos, niveles, ...
    'VariableNames', {'sesion','paciente','mes','condicion','completa','archivos','niveles'});

F = disenar_filtros(cfg);
pendientes = find(~cellfun(@(s) isfile(fullfile(dirS, [s '.mat'])), sesiones.sesion));
registrar(cfg, 'Etapa 1: %d sesiones (%d completas); %d por preprocesar.', ns, nnz(completa), numel(pendientes));
nw = usar_paralelo(cfg);
info = cell(numel(pendientes), 1);
arch = sesiones.archivos(pendientes);
nom  = sesiones.sesion(pendientes);
niv  = sesiones.niveles(pendientes);
parfor (i = 1:numel(pendientes), nw)
    t0 = tic;
    S = preprocesar_sesion(arch{i}, cfg, F);
    S.archivos = arch{i};
    S.niveles = niv{i};
    guardar_senal(fullfile(dirS, [nom{i} '.mat']), S);
    info{i} = S.info;
    fprintf('  %s: %.1f s, ICA iter=%d conv=%d, removidas=%s\n', nom{i}, toc(t0), ...
        S.info.iteraciones, S.info.convergio, mat2str(S.info.removidos));
end
for i = 1:numel(pendientes)
    registrar(cfg, '  %s: ICA iter=%d, convergió=%d, componentes removidas=%s, ASR=%s %%', nom{i}, ...
        info{i}.iteraciones, info{i}.convergio, mat2str(info{i}.removidos), ...
        mat2str(round(100 * info{i}.frac_asr, 1)));
end
sal = sesiones;
sal.archivos = cellfun(@(c) strjoin(cellfun(@nombre_corto, c, 'UniformOutput', false), ';'), ...
    sal.archivos, 'UniformOutput', false);
sal.niveles = cellfun(@mat2str, sal.niveles, 'UniformOutput', false);
writetable(sal, fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'Encoding', 'UTF-8');
end

function guardar_senal(ruta, S)
% Escritura atómica: si se interrumpe, no queda un .mat a medias.
tmp = [ruta '.tmp.mat'];
save(tmp, '-struct', 'S', '-v7');
movefile(tmp, ruta, 'f');
end

function n = nombre_corto(r)
[~, n] = fileparts(r);
end
