%% CORRIDA_CAR  Re-ejecuta el análisis con referencia promedio común (CAR) y
%  rechazo de artefactos, en resultados_car/ (los resultados previos sin CAR
%  quedan intactos en resultados/ para comparar).
%
%  Pasos: etapas 1-3 (512 muestras), etapa 4 fase A (E9, 512), verificación
%  de correlaciones entre canales tras la limpieza, etapa 7 (calibración),
%  etapa 8 (reposo/eventos) y los controles de sesgo (herencia del reposo,
%  control de mitades, orden invertido).
%  Los hiperparámetros de los clasificadores se reutilizan de la corrida sin
%  CAR (se copian), para que la comparación aísle el efecto del preprocesamiento.

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_car');
cfg = config_niveles();
if ~isfolder(cfg.rutas.resultados), mkdir(cfg.rutas.resultados); end
for sub = {'ajuste', 'ajuste_calibrado'}
    d = fullfile(cfg.rutas.resultados, sub{1});
    if ~isfolder(d), mkdir(d); end
    copyfile(fullfile(raiz, 'resultados', sub{1}, 'hiperparametros_elegidos.csv'), d);
end
cfg.segmentacion.duraciones = 512;
cfg.experimentos.fase = 'A';
cfg.calibracion.ajustar_despues = false;
save(fullfile(cfg.rutas.resultados, 'cfg_usada.mat'), 'cfg');
registrar(cfg, '==== Corrida con CAR + rechazo de artefactos (umbral %g MAD, margen %g s) ====', ...
    cfg.artefactos.umbral_z, cfg.artefactos.margen_s);

etapa1_preprocesar(cfg);

%% Verificación tras la limpieza: pares de canales casi idénticos y rechazo
dirS = fullfile(cfg.rutas.resultados, 'senales');
f = dir(fullfile(dirS, '*.mat'));
filas = {};
for k = 1:numel(f)
    S = load(fullfile(f(k).folder, f(k).name));
    for j = 1:numel(S.X)
        R = corrcoef(double(S.X{j}).'); R(1:size(R,1)+1:end) = 0;
        [~, a] = fileparts(S.archivos{j});
        filas(end+1, :) = {erase(f(k).name, '.mat'), a, nnz(triu(abs(R) > 0.98)), ...
            max(abs(R(:))), 100 * S.info.frac_mala(j)}; %#ok<AGROW>
    end
end
V = cell2table(filas, 'VariableNames', {'sesion','archivo','pares_r_mayor_098','r_max','pct_muestras_rechazadas'});
writetable(V, fullfile(cfg.rutas.resultados, 'verificacion_limpieza.csv'));
registrar(cfg, 'Tras CAR+ICA: pares r>0.98 por archivo = %.1f (máx %d); muestras rechazadas = %.1f %% (máx %.1f %%).', ...
    mean(V.pares_r_mayor_098), max(V.pares_r_mayor_098), mean(V.pct_muestras_rechazadas), max(V.pct_muestras_rechazadas));

etapa2_caracteristicas(cfg);
etapa3_exportar_arff(cfg);
etapa4_clasificar_weka(cfg);
etapa7_calibrar_caracterizadores(cfg);
etapa8_reposo_eventos(cfg);
analisis_confundido_reposo;
analisis_deriva_temporal;
prueba_orden_invertido;
registrar(cfg, '==== Corrida CAR terminada ====');
