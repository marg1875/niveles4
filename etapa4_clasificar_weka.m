function R = etapa4_clasificar_weka(cfg)
%ETAPA4_CLASIFICAR_WEKA  Clasificación en Weka con diseño experimental escalonado.
%   Fase A (duración): todas las duraciones x configuración de referencia
%     (cfg.canales.config_referencia, E9 = 16 canales) x rejilla completa.
%   Mejor duración: criterio PREFIJADO en cfg.experimentos.criterio
%     (máxima precisión media intra-sujeto en KMI 3c con HRS64+HO+HV,
%     promediada entre clasificadores; empate -> la duración más larga).
%   Fase B (canales): las 14 configuraciones x rejilla completa en la
%     duración de referencia (512) y en la mejor duración de la fase A.
%   Fuga: las cfg.weka.top_fuga mejores combinaciones intra-sujeto se
%     reevalúan con CV agrupada por bloque.
%   Rejilla completa = 10 conjuntos x {KMI,FIS} x {3c,2c} x {intra,mes,global}
%     x clasificadores.
%   Reanudable: cada dataset guarda su caché en resultados/weka/cache/.
%   -> resultados/weka/metricas.csv, resultados/weka/mejor_duracion.txt

dirW = fullfile(cfg.rutas.resultados, 'weka');
if ~isfolder(dirW), mkdir(dirW); end
clfs = {'kNN', 'MLP', 'SVM'};
if cfg.weka.incluir_extra, clfs = [clfs, {'BN','NB','RF','RT'}]; end
conjuntos = {definir_conjuntos(cfg).nombre};
configs = {cfg.canales.configs.id};
if cfg.modo_prueba
    configs = {cfg.canales.config_referencia, 'E6'};
end
duraciones = cfg.segmentacion.duraciones;
dref = cfg.segmentacion.duracion_referencia;

if contains(cfg.experimentos.fase, 'A')
    registrar(cfg, 'Etapa 4 - fase A: %d duraciones x config %s.', numel(duraciones), cfg.canales.config_referencia);
    ejecutar(cfg, duraciones, {cfg.canales.config_referencia}, conjuntos, clfs, 'estratificada', [], dirW);
end
R = consolidar(dirW);
mejor = mejor_duracion(R, cfg);
if ~isnan(mejor)
    registrar(cfg, 'Mejor duración (criterio prefijado): %d muestras = %.3f s.', mejor, mejor/cfg.fs);
    fid = fopen(fullfile(dirW, 'mejor_duracion.txt'), 'w'); fprintf(fid, '%d\n', mejor); fclose(fid);
end
if contains(cfg.experimentos.fase, 'B')
    dB = unique([dref, mejor(~isnan(mejor))]);
    registrar(cfg, 'Etapa 4 - fase B: %d configuraciones x duraciones %s.', numel(configs), mat2str(dB));
    ejecutar(cfg, dB, configs, conjuntos, clfs, 'estratificada', [], dirW);
    % Chequeo de fuga: mejores combinaciones intra (KMI 3c) con CV agrupada
    R = consolidar(dirW);
    top = mejores_combinaciones(R, cfg);
    if ~isempty(top)
        registrar(cfg, 'Chequeo de fuga (CV agrupada) en %d combinaciones.', height(top));
        for i = 1:height(top)
            ejecutar(cfg, top.duracion(i), top.config(i), top.conjunto(i), top.clasificador(i), ...
                'agrupada', {'intra'}, dirW);
        end
    end
end
R = consolidar(dirW);
writetable(R, fullfile(dirW, 'metricas.csv'), 'Encoding', 'UTF-8');
registrar(cfg, 'Etapa 4: %d resultados (%d ok) en weka/metricas.csv.', height(R), nnz(strcmp(R.estado,'ok')));
end

function ejecutar(cfg, duraciones, configs, conjuntos, clfs, cv_modo, alcances, dirW)
% Construye las unidades de trabajo (una por duración x dataset) y las
% ejecuta en paralelo; cada unidad evalúa todas sus combinaciones.
[cc, cj, cl] = ndgrid(1:numel(configs), 1:numel(conjuntos), 1:numel(clfs));
col = @(x, i) reshape(x(i(:)), [], 1);
combos = table(col(configs, cc), col(conjuntos, cj), col(clfs, cl), ...
    repmat({cv_modo}, numel(cc), 1), 'VariableNames', {'config','conjunto','clasificador','cv_modo'});
unidades = struct('arff', {}, 'cache', {}, 'info', {});
for L = duraciones(:).'
    dirA = fullfile(cfg.rutas.resultados, 'arff', sprintf('d%d', L));
    if ~isfile(fullfile(dirA, 'indice.csv')), continue, end
    I = readtable(fullfile(dirA, 'indice.csv'), 'TextType', 'char');
    if ~isempty(alcances), I = I(ismember(I.alcance, alcances), :); end
    dirCache = fullfile(dirW, 'cache', sprintf('d%d', L));
    if ~isfolder(dirCache), mkdir(dirCache); end
    for k = 1:height(I)
        u.arff  = fullfile(dirA, I.archivo{k});
        u.cache = fullfile(dirCache, [I.dataset{k} '.mat']);
        u.info  = struct('duracion', L, 'dataset', I.dataset{k}, 'condicion', I.condicion{k}, ...
            'esquema', I.esquema{k}, 'alcance', I.alcance{k}, 'id', I.id{k});
        unidades(end+1) = u; %#ok<AGROW>
    end
end
% Los datasets globales (los más lentos) primero, para equilibrar la carga
orden = [find(strcmp(arrayfun(@(u) u.info.alcance, unidades, 'UniformOutput', false), 'global')), ...
         find(strcmp(arrayfun(@(u) u.info.alcance, unidades, 'UniformOutput', false), 'mes')), ...
         find(strcmp(arrayfun(@(u) u.info.alcance, unidades, 'UniformOutput', false), 'intra'))];
unidades = unidades(orden);
registrar(cfg, '  %d datasets x %d combinaciones (CV %s).', numel(unidades), height(combos), cv_modo);
nw = usar_paralelo(cfg);
t0 = tic;
parfor (i = 1:numel(unidades), nw)
    evaluar_dataset(unidades(i).arff, combos, cfg, unidades(i).cache, unidades(i).info);
end
registrar(cfg, '  terminado en %.1f min.', toc(t0)/60);
end

function R = consolidar(dirW)
f = dir(fullfile(dirW, 'cache', '**', '*.mat'));
f = f(~endsWith({f.name}, '.tmp.mat'));
partes = cell(numel(f), 1);
for k = 1:numel(f)
    S = load(fullfile(f(k).folder, f(k).name), 'R');
    partes{k} = S.R;
end
R = vertcat(partes{:});
end

function mejor = mejor_duracion(R, cfg)
mejor = NaN;
if isempty(R), return, end
c = cfg.experimentos.criterio;
sel = strcmp(R.condicion, c.condicion) & strcmp(R.esquema, c.esquema) & ...
      strcmp(R.alcance, c.alcance) & strcmp(R.conjunto, c.conjunto) & ...
      strcmp(R.config, cfg.canales.config_referencia) & strcmp(R.cv_modo, 'estratificada') & ...
      strcmp(R.estado, 'ok');
if ~any(sel), return, end
[g, d] = findgroups(R.duracion(sel));
m = splitapply(@mean, R.acc(sel), g);
cand = d(abs(m - max(m)) < 1e-9);
mejor = max(cand);
end

function top = mejores_combinaciones(R, cfg)
c = cfg.experimentos.criterio;
sel = strcmp(R.condicion, c.condicion) & strcmp(R.esquema, c.esquema) & ...
      strcmp(R.alcance, 'intra') & strcmp(R.cv_modo, 'estratificada') & strcmp(R.estado, 'ok');
if ~any(sel), top = table(); return, end
[g, d, cf, cj, cl] = findgroups(R.duracion(sel), R.config(sel), R.conjunto(sel), R.clasificador(sel));
m = splitapply(@mean, R.acc(sel), g);
top = table(d, cf, cj, cl, m, 'VariableNames', {'duracion','config','conjunto','clasificador','acc_media'});
top = sortrows(top, 'acc_media', 'descend');
top = top(1:min(cfg.weka.top_fuga, height(top)), :);
end
