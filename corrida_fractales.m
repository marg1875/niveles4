%% CORRIDA_FRACTALES  Etapa 9 sobre el preprocesamiento final (resultados_asr/).
%  ASR k = 30 + CAR + ICA, sin rechazo; base = par calibrado de la etapa 7.
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
etapa9_fractales(cfg);
