function [M, meta] = cache_ventanas_kmi(arch)
%CACHE_VENTANAS_KMI  Ventanas contiguas KMI de un archivo w_*.mat, cargadas una vez por proceso.
%   [M, meta] = cache_ventanas_kmi(arch): M (ventanas x ancho, double) con
%   las curvas/valores de caract_curvas y meta (tabla sin 'feat'). Cada
%   worker del pool lee el archivo la primera vez y lo conserva en memoria
%   (variable persistente), en lugar de recibir la matriz serializada
%   desde el cliente (parallel.pool.Constant falla con matrices de cientos de MB).
persistent ruta Mc metac
if isempty(ruta) || ~strcmp(ruta, arch)
    S = load(arch);
    k = strcmp(S.W.condicion, 'KMI');
    Mc = double(S.W.feat(k, :));
    metac = removevars(S.W(k, :), 'feat');
    ruta = arch;
end
M = Mc; meta = metac;
end
