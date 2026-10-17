function m = especificar_modelo(clf, params, sel)
%ESPECIFICAR_MODELO  Especificación de modelo con identificador único.
%   m = especificar_modelo(clf, params, sel) devuelve un struct con clf,
%   params, sel, miembros (vacío) e id, p. ej.
%   'SVM|C=10;gamma=0.01;kernel=rbf|ig16'. El id se usa como clave de caché.
if nargin < 3, sel = 'ninguna'; end
m.clf = clf;
m.params = params;
m.sel = sel;
m.miembros = [];
m.id = sprintf('%s|%s|%s', clf, texto_params(params), sel);
end

function t = texto_params(p)
f = sort(fieldnames(p));
partes = cell(1, numel(f));
for i = 1:numel(f)
    v = p.(f{i});
    if isnumeric(v) || islogical(v), v = num2str(v); end
    partes{i} = sprintf('%s=%s', f{i}, v);
end
t = strjoin(partes, ';');
end
