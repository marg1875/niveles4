%% EJECUTAR_PIPELINE  Pipeline completo niveles4 (MATLAB + Weka), por lotes.
%
%   Reproduce la metodología de Martinez-Peon et al. (2024, J. Neural Eng.)
%   sobre los registros g.Nautilus de pacientes con ECV (carpeta Datos/):
%     1) Inventario + preprocesamiento por sesión (filtro, z-score, blanqueo,
%        ICA de gradiente natural con una W por sesión)  -> resultados/senales
%     2) Segmentación multi-duración + Hurst HO/HRS/HV  -> resultados/caracteristicas
%     3) Exportación ARFF para Weka                      -> resultados/arff
%     4) Clasificación en Weka (kNN, MLP, SVM; CV 5-fold) -> resultados/weka
%     5) Reportes (tablas y figuras)                     -> resultados/reportes
%     6) Ajuste de clasificadores, KMI 3 clases          -> resultados/ajuste
%
%   Todo es reanudable: cada etapa omite lo ya calculado. Para recalcular
%   una etapa, borre su carpeta en resultados/.
%   Uso:  >> ejecutar_pipeline                (corrida completa)
%         Modo de prueba: poner cfg.modo_prueba = true abajo.

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));

cfg = config_niveles();
% --- Ajustes locales de la corrida (opcional) ---------------------------
% cfg.modo_prueba = true;
% cfg.experimentos.fase = 'A';
% cfg.weka.incluir_extra = true;

if cfg.modo_prueba
    cfg.rutas.resultados = fullfile(raiz, 'resultados_prueba');
    cfg.segmentacion.duraciones = [512 125];
end
if ~isfolder(cfg.rutas.resultados), mkdir(cfg.rutas.resultados); end
save(fullfile(cfg.rutas.resultados, 'cfg_usada.mat'), 'cfg');
registrar(cfg, '==== niveles4 %s | modo_prueba=%d | etapas %s ====', cfg.version, ...
    cfg.modo_prueba, mat2str(cfg.etapas));

t0 = tic;
if ismember(1, cfg.etapas), etapa1_preprocesar(cfg); end
if ismember(2, cfg.etapas), etapa2_caracteristicas(cfg); end
if ismember(3, cfg.etapas), etapa3_exportar_arff(cfg); end
if ismember(4, cfg.etapas), etapa4_clasificar_weka(cfg); end
if ismember(5, cfg.etapas), etapa5_reportes(cfg); end
if ismember(6, cfg.etapas), etapa6_ajuste_clasificadores(cfg); end
registrar(cfg, '==== Pipeline terminado en %.1f min ====', toc(t0)/60);
