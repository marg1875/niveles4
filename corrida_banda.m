%% CORRIDA_BANDA  Etapa 16 (potencia de banda optimizada y combinada) con los 56 registros. Reanudable.
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
etapa16_banda_combinada(cfg);
