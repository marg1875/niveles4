%% CORRIDA_REVISION  Etapa 15 (controles de revisión) sobre resultados_asr/. Reanudable.
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
etapa15_revision(cfg);
