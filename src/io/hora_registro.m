function h = hora_registro(ruta)
%HORA_REGISTRO  Fecha y hora de fin de grabación de un .mat ("Created on" del encabezado).
%   h = hora_registro(ruta) devuelve un datetime. El encabezado de texto de
%   los .mat v5 trae la marca 'Created on: Mon Oct 12 10:31:05 2026', que
%   corresponde al momento en que se guardó el registro.
fid = fopen(ruta, 'r');
if fid < 0, error('hora_registro:abrir', 'No se pudo abrir %s', ruta); end
b = fread(fid, 116, '*char')';
fclose(fid);
tok = regexp(b, 'Created on: (\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})', 'tokens', 'once');
if isempty(tok), h = NaT; return, end
h = datetime(regexprep(tok{1}, ' +', ' '), 'InputFormat', 'eee MMM d HH:mm:ss yyyy', 'Locale', 'en_US');
end
