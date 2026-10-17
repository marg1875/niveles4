%% CORRIDA_NIVELES3_GANADORES  Etapa 11 sobre resultados_asr/ (requiere etapas 9 y 10).
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
etapa11_niveles3_ganadores(cfg);
