%% CORRIDA_HIG_DFA  Etapa 10 (calibración de Higuchi y DFA) sobre resultados_asr/.
raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
etapa10_calibrar_hig_dfa(cfg);
