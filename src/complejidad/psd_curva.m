function [P, f] = psd_curva(V, fs, nfft)
%PSD_CURVA  Densidad espectral de Welch para el exponente aperiódico 1/f.
%   [P, f] = psd_curva(V, fs, nfft), V (C x N) -> P (C x nf) en f (Hz).
%   Ventana de Hann de nfft muestras con 50 % de traslape (si N < nfft, una
%   sola ventana de N). La pendiente de log P vs log f en una banda da
%   -beta (P ~ 1/f^beta); beta se relaciona con el exponente de escala
%   (beta = 2H + 1 para fBm). No es una potencia de banda: mide la forma
%   del fondo 1/f, independiente de la amplitud.
N = size(V, 2);
w = min(nfft, N);
[P, f] = pwelch(double(V).', hann(w), floor(w/2), max(nfft, w), fs);
P = P.'; f = f.';
P(P <= 0) = NaN;
end
