function D = cargar_datos_ajuste(ruta, cfg, conjunto)
%CARGAR_DATOS_AJUSTE  Carga un ARFF y conserva la configuración de canales y el conjunto.
%   D = cargar_datos_ajuste(ruta, cfg) usa cfg.ajuste.config y
%   cfg.ajuste.conjunto; D = cargar_datos_ajuste(ruta, cfg, conjunto) usa
%   otro conjunto de características (nombre de definir_conjuntos o texto
%   'HRS64+HO', etc.: familias separadas por '+').
if nargin < 3, conjunto = cfg.ajuste.conjunto; end
weka_iniciar(cfg);
data = weka_cargar_arff(ruta);
canales = cfg.canales.configs(strcmp({cfg.canales.configs.id}, cfg.ajuste.config)).canales;
metodos = strsplit(conjunto, '+');
attrs = {};
for m = 1:numel(metodos), attrs = [attrs, strcat(metodos{m}, '_', canales)]; end %#ok<AGROW>
D = weka_seleccionar_atributos(data, attrs);
end
