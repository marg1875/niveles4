function ms = niveles_rangors(N, p, bloque_min)
%NIVELES_RANGORS  Tamaños de bloque que usa rangors(v, p) para una serie de N.
%   ms = niveles_rangors(N, p) reproduce la regla de rangors.m:
%       while m >= N/p,  m = floor(N/2^i);  ...;  i = i+1;  end
%   es decir, m_i = floor(N/2^i) para i = 0,1,... hasta el PRIMER m_i < N/p
%   (incluido). Ej.: N = 512, p = 64 -> m = 512,256,...,8,4 (8 niveles).
%   Se descartan los niveles con bloques menores que bloque_min (por
%   defecto 2: un bloque de 1 muestra tiene S = 0).
if nargin < 3, bloque_min = 2; end
ms = [];
m = N; i = 0;
while m >= N/p
    m = floor(N / 2^i);
    ms(end+1) = m; %#ok<AGROW>
    i = i + 1;
    if m == 0, break, end
end
ms = ms(ms >= bloque_min);
end
