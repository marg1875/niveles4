function R = buscar_en_sesion(ruta, modelos, cfg, ruta_cache)
%BUSCAR_EN_SESION  Evalúa (CV 5-fold, semilla cfg.weka.semilla) cada modelo en una sesión.
%   R = buscar_en_sesion(ruta, modelos, cfg, ruta_cache) devuelve una tabla
%   con id, clf, sel, acc, kappa, mae, wf y estado. Reanudable: los modelos
%   ya presentes en ruta_cache no se recalculan, y la caché se guarda cada
%   25 modelos (escritura atómica). Pensada para un parfor por sesión.
R = table();
if isfile(ruta_cache)
    S = load(ruta_cache, 'R'); R = S.R;
end
hechos = {};
if ~isempty(R), hechos = R.id; end
pend = modelos(~ismember({modelos.id}, hechos));
if isempty(pend), return, end
D = cargar_datos_ajuste(ruta, cfg);
nuevas = cell(numel(pend), 1);
ultimo = 0;                                   % último índice ya volcado a R
for i = 1:numel(pend)
    m = pend(i);
    fila = table({m.id}, {m.clf}, {m.sel}, NaN, NaN, NaN, NaN, {'ok'}, ...
        'VariableNames', {'id','clf','sel','acc','kappa','mae','wf','estado'});
    try
        M = weka_validar(D, weka_clasificador_param(m, cfg), cfg);
        fila.acc = M.acc; fila.kappa = M.kappa; fila.mae = M.mae; fila.wf = M.wf;
    catch err
        fila.estado = {['error: ' err.message]};
    end
    nuevas{i} = fila;
    if mod(i, 25) == 0 || i == numel(pend)
        R = [R; vertcat(nuevas{ultimo+1:i})]; %#ok<AGROW>
        ultimo = i;
        guardar(ruta_cache, R);
        % avance legible por progreso.ps1: "<hechos> <total>"
        fid = fopen(strrep(ruta_cache, '.mat', '.avance.txt'), 'w');
        if fid > 0, fprintf(fid, '%d %d\n', height(R), height(R) + numel(pend) - i); fclose(fid); end
    end
end
end

function guardar(ruta, R)
tmp = [ruta '.tmp.mat'];
save(tmp, 'R');
movefile(tmp, ruta, 'f');
end
