function nombres = nombres_caracteristicas(cfg)
%NOMBRES_CARACTERISTICAS  Nombres de atributos en el orden de extraer_hurst_ventana.
%   'HO_<canal>' ..., 'HRS<p>_<canal>' ... (para cada p), 'HV_<canal>' ...
ch = cfg.canales.nombres;
metodos = [{'HO'}, arrayfun(@(p) sprintf('HRS%d', p), cfg.hurst.particiones_hrs, 'UniformOutput', false), {'HV'}];
nombres = {};
for m = 1:numel(metodos)
    nombres = [nombres, strcat(metodos{m}, '_', ch)]; %#ok<AGROW>
end
end
