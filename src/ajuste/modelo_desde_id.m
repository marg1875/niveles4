function m = modelo_desde_id(id, sel)
%MODELO_DESDE_ID  Reconstruye una especificación de modelo a partir de su id.
%   m = modelo_desde_id('SVM|C=1;E=1;filtro=0;kernel=poly|ig24') invierte
%   especificar_modelo. Con un segundo argumento se reemplaza la selección
%   de características (p. ej. 'ninguna').
p = strsplit(id, '|');
params = struct();
if numel(p) >= 2 && ~isempty(p{2})
    for kv = strsplit(p{2}, ';')
        par = strsplit(kv{1}, '=');
        v = str2double(par{2});
        if isnan(v) || any(par{2} == ',') || strcmp(par{1}, 'H'), v = par{2}; end   % texto
        params.(par{1}) = v;
    end
end
if nargin < 2, sel = p{3}; end
m = especificar_modelo(p{1}, params, sel);
end
