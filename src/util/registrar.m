function registrar(cfg, formato, varargin)
%REGISTRAR  Imprime un mensaje con marca de tiempo y lo añade al log.
%   registrar(cfg, formato, ...) escribe en consola y en
%   resultados/pipeline.log (mismo uso que fprintf).
msg = sprintf(formato, varargin{:});
linea = sprintf('[%s] %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg); %#ok<TNOW1,DATST>
fprintf('%s\n', linea);
if ~isfolder(cfg.rutas.resultados), mkdir(cfg.rutas.resultados); end
fid = fopen(fullfile(cfg.rutas.resultados, 'pipeline.log'), 'a', 'n', 'UTF-8');
if fid > 0
    fprintf(fid, '%s\n', linea);
    fclose(fid);
end
end
