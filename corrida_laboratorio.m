%% CORRIDA_LABORATORIO  Etapa 12 (laboratorio Higuchi/DFA) sobre resultados_asr/.
%  Reanudable: si se interrumpe, al volver a correrla retoma desde la caché.
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
etapa12_laboratorio_hig_dfa(cfg);
