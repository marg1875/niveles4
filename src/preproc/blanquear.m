function [Z, Q, info] = blanquear(X, cfg)
%BLANQUEAR  Blanqueo PCA por eigendescomposición (Ecs. 2-3 del artículo).
%   [Z, Q] = blanquear(X, cfg) con X (m canales x n muestras, media cero).
%   Q (k x m) es la matriz de blanqueo: Z = Q*X tiene cov(Z) = I (modo
%   'estandar'). Se descartan los autovalores < cfg.blanqueo.tol_rango*max
%   (deficiencia de rango, p. ej. por referencia común), por lo que k <= m.
%
%   Modo 'p0_articulo': aplica además el operador diagonal P0 de la Ec. 3
%   (pesos para los 2 autovalores mayores, los intermedios y los 2 menores,
%   cfg.blanqueo.p0_pesos). El artículo da rangos (0.7-0.9 / 0.3-0.7 /
%   0.1-0.3); por defecto se usan los puntos medios.
X = X - mean(X, 2);
C = (X * X.') / size(X, 2);
C = (C + C.')/2;
[E, D] = eig(C);
[lam, orden] = sort(diag(D), 'descend');
E = E(:, orden);
mantener = lam > cfg.blanqueo.tol_rango * lam(1);
lam = lam(mantener);
E = E(:, mantener);
Q = diag(1./sqrt(lam)) * E.';

if strcmp(cfg.blanqueo.modo, 'p0_articulo')
    k = numel(lam);
    p0 = cfg.blanqueo.p0_pesos(2) * ones(k,1);
    p0(1:min(2,k)) = cfg.blanqueo.p0_pesos(1);
    if k > 4, p0(end-1:end) = cfg.blanqueo.p0_pesos(3); end
    Q = diag(p0) * Q;
end
Z = Q * X;
info.autovalores = lam;
info.rango = numel(lam);
end
