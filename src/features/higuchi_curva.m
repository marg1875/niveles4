function Lk = higuchi_curva(V, kmax)
%HIGUCHI_CURVA  Longitud media de curva de Higuchi L(k), k = 1..kmax.
%   Lk = higuchi_curva(V, kmax), V (C x N) -> Lk (C x kmax). Para cada
%   inicio m = 1..k (base 1): L_m(k) = sum |x(m+ik) - x(m+(i-1)k)| *
%   (N-1) / (floor((N-m+1)/k) k), omitiendo subseries con <= 2 puntos;
%   L(k) = mean_m L_m(k) / k (normalización de Niveles3/compute_higuchi).
%   Valores <= 1e-10 o no calculables -> NaN.
[C, N] = size(V);
X = double(V);
Lk = nan(C, kmax);
for k = 1:kmax
    suma = zeros(C, 1); validos = 0;
    for m = 1:k
        idx = m:k:N;
        if numel(idx) <= 2, continue, end
        largo = sum(abs(diff(X(:, idx), 1, 2)), 2);
        suma = suma + largo * (N - 1) / (floor((N - m + 1) / k) * k);
        validos = validos + 1;
    end
    if validos > 0, Lk(:, k) = suma / (k * validos); end
end
Lk(Lk <= 1e-10) = NaN;
end
