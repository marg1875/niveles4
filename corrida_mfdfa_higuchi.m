%% CORRIDA_MFDFA_HIGUCHI  Etapa 14 sobre resultados_asr/ (requiere la etapa 13). Reanudable.
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
etapa14_mfdfa_higuchi(cfg);
