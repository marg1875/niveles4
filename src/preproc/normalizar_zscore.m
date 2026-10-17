function Y = normalizar_zscore(X)
%NORMALIZAR_ZSCORE  Normalización por canal (Ec. 1 del artículo):
%   X_i(t) = (x_i(t) - mu_i) / sigma_i, con X de tamaño canales x muestras.
mu = mean(X, 2);
sd = std(X, 0, 2);
sd(sd < eps) = 1;
Y = (X - mu) ./ sd;
end
