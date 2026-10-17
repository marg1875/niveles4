function [obs, nula, p] = permutacion_secuencial(T, nombres, modelo, cfg, nw, nmax, h, archAvance)
%PERMUTACION_SECUENCIAL  Prueba de permutación de la LOSO 10 vs 40 con paro temprano.
%   [obs, nula, p] = permutacion_secuencial(T, nombres, modelo, cfg, nw, nmax, h, archAvance)
%   Igual que la permutación por sesión completa usada hasta ahora
%   (etiquetas 10/40 intercambiadas en sesiones al azar, semillas 500+i),
%   pero con el procedimiento secuencial de Besag y Clifford (1991): se
%   detiene en cuanto h permutaciones igualan o superan el valor observado
%   (p = h / n), o al llegar a nmax (p = (1 + excedentes) / (1 + nmax)).
%   Un resultado claramente no significativo termina con ~h/p permutaciones
%   en lugar de nmax; los significativos se calculan completos.
%   archAvance (opcional): archivo de avance 'hechas nmax'.
if nargin < 7 || isempty(h), h = 20; end
if nargin < 8, archAvance = ''; end
obs = mean(loso_niveles(T, nombres, modelo, cfg, []));
sesC = unique(T.sesion);
nula = nan(nmax, 1);
lote = max(1, nw) * 2;
n = 0;
while n < nmax
    idx = n+1 : min(nmax, n + lote);
    parte = zeros(numel(idx), 1);
    parfor (j = 1:numel(idx), nw)
        rs = RandStream('mt19937ar', 'Seed', 500 + idx(j));
        parte(j) = mean(loso_niveles(T, nombres, modelo, cfg, sesC(rand(rs, numel(sesC), 1) < 0.5)));
    end
    nula(idx) = parte;
    n = idx(end);
    if ~isempty(archAvance)
        fid = fopen(archAvance, 'w'); if fid > 0, fprintf(fid, '%d %d', n, nmax); fclose(fid); end
    end
    if nnz(nula(1:n) >= obs) >= h, break, end
end
nula = nula(1:n);
e = nnz(nula >= obs);
if e >= h, p = e / n; else, p = (1 + e) / (1 + n); end
end
