function s = ventanas_reposo_inicial(L, cfg)
%VENTANAS_REPOSO_INICIAL  Inicios de ventanas no traslapadas de L muestras
%   dentro del reposo inicial [cfg.eventos.reposo_desde_s, cfg.eventos.reposo_inicial_s].
fs = cfg.fs;
a = round(cfg.eventos.reposo_desde_s * fs) + 1;
b = round(cfg.eventos.reposo_inicial_s * fs);
n = floor((b - a + 1) / L);
s = a + (0:n-1) * L;
end
