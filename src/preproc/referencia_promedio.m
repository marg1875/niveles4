function Y = referencia_promedio(X)
%REFERENCIA_PROMEDIO  Referencia promedio común (CAR).
%   Y = referencia_promedio(X), X (canales x muestras): a cada muestra se le
%   resta la media de todos los canales. Elimina la señal común a todos los
%   electrodos (contaminación de la referencia, artefactos globales). El
%   resultado tiene rango (canales - 1), lo que el blanqueo ya contempla.
Y = X - mean(X, 1, 'omitnan');               % los canales eliminados (NaN) no cuentan
end
