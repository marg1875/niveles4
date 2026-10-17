%% CORRIDA_ASR  Análisis principal con el preprocesamiento final:
%  filtro -> ASR (k = 30) -> CAR -> ICA por sesión, SIN rechazo de ventanas
%  y con todas las sesiones, en resultados_asr/.
%
%  Las señales son las de comparar_asr.m (resultados_var/asr_k30/senales),
%  con la misma configuración; se copian quitando las máscaras de rechazo
%  (aquella corrida aún las marcaba) para no repetir el ICA.
%  Pasos: etapas 1-4 (512 muestras, fase A), verificación de limpieza,
%  etapa 7 (calibración), etapa 8 (reposo/eventos) y los controles de sesgo
%  (herencia del reposo, control de mitades, orden invertido).
%  Los hiperparámetros de los clasificadores se reutilizan de resultados/.

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
assert(cfg.asr.activar && cfg.asr.k == 30 && ~cfg.artefactos.activar, 'config_niveles: se espera ASR k = 30 sin rechazo');
dirS = fullfile(cfg.rutas.resultados, 'senales');
if ~isfolder(dirS), mkdir(dirS); end
for sub = {'ajuste', 'ajuste_calibrado'}
    d = fullfile(cfg.rutas.resultados, sub{1});
    if ~isfolder(d), mkdir(d); end
    copyfile(fullfile(raiz, 'resultados', sub{1}, 'hiperparametros_elegidos.csv'), d);
end
cfg.segmentacion.duraciones = 512;
cfg.experimentos.fase = 'A';
cfg.calibracion.ajustar_despues = false;
save(fullfile(cfg.rutas.resultados, 'cfg_usada.mat'), 'cfg');
registrar(cfg, '==== Corrida final: ASR k = %d + CAR + ICA, sin rechazo de ventanas ====', cfg.asr.k);

%% Señales con ASR k = 30 (sin máscaras de rechazo)
origen = fullfile(raiz, 'resultados_var', 'asr_k30', 'senales');
f = dir(fullfile(origen, '*.mat'));
for k = 1:numel(f)
    destino = fullfile(dirS, f(k).name);
    if isfile(destino), continue, end
    S = load(fullfile(f(k).folder, f(k).name));
    S.malas = repmat({zeros(0, 2)}, 1, numel(S.X));
    S.info.frac_mala = zeros(1, numel(S.X));
    save(destino, '-struct', 'S', '-v7');
end
etapa1_preprocesar(cfg);                           % inventario y sesiones (las señales ya existen)

%% Verificación tras la limpieza: pares de canales casi idénticos y % reparado
filas = {};
for k = 1:numel(f)
    S = load(fullfile(dirS, f(k).name));
    for j = 1:numel(S.X)
        R = corrcoef(double(S.X{j}).'); R(1:size(R,1)+1:end) = 0;
        [~, a] = fileparts(S.archivos{j});
        filas(end+1, :) = {erase(f(k).name, '.mat'), a, nnz(triu(abs(R) > 0.98)), ...
            max(abs(R(:))), 100 * S.info.frac_asr(j)}; %#ok<AGROW>
    end
end
V = cell2table(filas, 'VariableNames', {'sesion','archivo','pares_r_mayor_098','r_max','pct_reparado_asr'});
writetable(V, fullfile(cfg.rutas.resultados, 'verificacion_limpieza.csv'));
registrar(cfg, 'Tras ASR+CAR+ICA: pares r>0.98 por archivo = %.1f (máx %d); registro reparado por ASR = %.1f %% (máx %.1f %%).', ...
    mean(V.pares_r_mayor_098), max(V.pares_r_mayor_098), mean(V.pct_reparado_asr), max(V.pct_reparado_asr));

etapa2_caracteristicas(cfg);
etapa3_exportar_arff(cfg);
etapa4_clasificar_weka(cfg);
etapa7_calibrar_caracterizadores(cfg);
etapa8_reposo_eventos(cfg);
analisis_confundido_reposo;
analisis_deriva_temporal;
prueba_orden_invertido;
registrar(cfg, '==== Corrida final terminada ====');
