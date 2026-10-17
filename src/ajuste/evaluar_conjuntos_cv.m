function evaluar_conjuntos_cv(ruta, conjuntos, modelos, etiq, cfg, ruta_cache)
%EVALUAR_CONJUNTOS_CV  CV de 5 folds de varios conjuntos de familias, con caché.
%   evaluar_conjuntos_cv(ruta, conjuntos, modelos, etiq, cfg, ruta_cache)
%   ruta: ARFF con atributos '<familia>_<canal>'. conjuntos: celda de
%   nombres 'FamA+FamB' (se usan todos los canales de cada familia).
%   modelos/etiq: especificaciones (modelo_desde_id) y sus etiquetas.
%   Añade a ruta_cache (tabla R: conjunto, clasificador, acc, kappa) los
%   pares conjunto x clasificador que falten (reanudable).
R = table();
if isfile(ruta_cache), S = load(ruta_cache, 'R'); R = S.R; end
hechos = {};
if ~isempty(R), hechos = strcat(R.conjunto, '|', R.clasificador); end
weka_iniciar(cfg);
data = weka_cargar_arff(ruta);
filas = {};
for c = 1:numel(conjuntos)
    if all(ismember(strcat(conjuntos{c}, '|', etiq), hechos)), continue, end
    fam = strsplit(conjuntos{c}, '+');
    attrs = {};
    for f = 1:numel(fam), attrs = [attrs, strcat(fam{f}, '_', cfg.canales.nombres)]; end %#ok<AGROW>
    D = weka_seleccionar_atributos(data, attrs);
    for i = 1:numel(modelos)
        M = weka_validar(D, weka_clasificador_param(modelos(i), cfg), cfg);
        filas(end+1, :) = {conjuntos{c}, etiq{i}, M.acc, M.kappa}; %#ok<AGROW>
    end
end
if isempty(filas), return, end
R = [R; cell2table(filas, 'VariableNames', {'conjunto','clasificador','acc','kappa'})];
tmp = [ruta_cache '.tmp.mat']; save(tmp, 'R'); movefile(tmp, ruta_cache, 'f');
end
