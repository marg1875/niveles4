function n = usar_paralelo(cfg)
%USAR_PARALELO  Abre (si procede) un pool paralelo y devuelve su tamaño.
%   n = usar_paralelo(cfg) devuelve 0 si no hay Parallel Computing Toolbox,
%   si cfg.paralelo.activar es false o si el pool no puede iniciarse; en ese
%   caso los parfor del pipeline se ejecutan en serie.
n = 0;
if ~cfg.paralelo.activar, return, end
try
    if ~license('test', 'Distrib_Computing_Toolbox') || isempty(ver('parallel'))
        return
    end
    p = gcp('nocreate');
    if isempty(p)
        % El perfil 'Processes' limita los workers a los núcleos físicos; se
        % amplía en memoria (sin guardar el perfil) para usar también los
        % hilos SMT (cfg.paralelo.workers).
        c = parcluster('Processes');
        if c.NumWorkers < cfg.paralelo.workers, c.NumWorkers = cfg.paralelo.workers; end
        p = parpool(c, cfg.paralelo.workers);
    end
    n = p.NumWorkers;
catch err
    warning('usar_paralelo:pool', 'Sin pool paralelo (%s). Se ejecuta en serie.', err.message);
    n = 0;
end
end
