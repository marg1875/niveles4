function T = cargar_cache_conjuntos(dirC, prefijo, ids)
%CARGAR_CACHE_CONJUNTOS  Une las cachés de evaluar_conjuntos_cv (una por sesión).
partes = cell(numel(ids), 1);
for s = 1:numel(ids)
    S = load(fullfile(dirC, [prefijo ids{s} '.mat']), 'R');
    R = S.R; R.sesion = repmat(ids(s), height(R), 1);
    partes{s} = R;
end
T = vertcat(partes{:});
end
