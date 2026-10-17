function G = semivariograma(V, hmax)
%SEMIVARIOGRAMA  Semivariograma empírico (Ec. 16 del artículo).
%   G = semivariograma(V, hmax), con V (C x N). G (C x hmax) con
%       gamma(h) = 1/(2 n(h)) * sum_{t=1}^{n(h)} (x(t+h) - x(t))^2,  n(h) = N-h
%   para h = 1..hmax (por defecto floor(N/2)).
%
%   Equivale al bucle de variograma.m del grupo, salvo que aquí se incluye
%   el lag 1 (el original empieza en i = 2 y deja vario(1) = 0).
%   Implementación vectorizada: sum (x_{t+h}-x_t)^2 = S2(t<=N-h) +
%   S2(t>h) - 2 r(h), con r(h) la autocorrelación (FFT) sin normalizar.
[C, N] = size(V);
if nargin < 2, hmax = floor(N/2); end
V = double(V);
h = 1:hmax;
s2 = cumsum(V.^2, 2);
nfft = 2^nextpow2(2*N);
F = fft(V, nfft, 2);
r = real(ifft(abs(F).^2, [], 2));
r = r(:, h + 1);                                     % r(h) = sum x_t x_{t+h}
suma = s2(:, N - h) + (s2(:, N) - s2(:, h)) - 2*r;
G = suma ./ (2 * (N - h));
G = max(G, 0);                                       % ruido de redondeo
if C == 0, G = zeros(0, hmax); end
end
