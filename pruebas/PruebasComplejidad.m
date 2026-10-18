classdef PruebasComplejidad < matlab.unittest.TestCase
%PRUEBASCOMPLEJIDAD  Algoritmos de complejidad/escalamiento de la etapa 13.
%   Referencias teóricas: ruido blanco (H = 0.5, PSD plana, PE ~ 1, LZC ~ 1),
%   movimiento browniano (alfa = 1.5, PSD ~ 1/f^2) y senoide (regular).
    properties
        R    % ruido blanco 8 x 2048
        B    % browniano
        S    % senoides
    end

    methods (TestClassSetup)
        function datos(tc)
            raiz = fileparts(fileparts(mfilename('fullpath')));
            addpath(raiz, genpath(fullfile(raiz, 'src')));
            rng(31);
            tc.R = randn(8, 2048);
            tc.B = cumsum(randn(8, 2048), 2);
            tc.S = sin(2*pi*(1:8).' * (0:2047) / 256);
        end
    end

    methods (Test)
        function mfdfaMonofractalYDfa(tc)
            ns = [16 24 32 48 64 96 128];
            F = mfdfa_curva(tc.R, ns, 1, [-3 2 3]);
            pend = @(f) arrayfun(@(c) polyfit(log10(ns), log10(f(c, :)), 1) * [1; 0], 1:8);
            h = [pend(F(:, :, 1)); pend(F(:, :, 2)); pend(F(:, :, 3))];
            tc.verifyEqual(mean(h, 2), 0.5 * ones(3, 1), 'AbsTol', 0.12);
            tc.verifyLessThan(mean(h(1, :) - h(3, :)), 0.2);         % ancho multifractal pequeño
            tc.verifyEqual(F(:, :, 2), dfa_curva(tc.R, ns, 1, 1), 'RelTol', 1e-10);   % q = 2 es el DFA
        end

        function dmaYAgregados(tc)
            ns = [9 15 25 41 65 101];
            p = @(F) mean(arrayfun(@(c) polyfit(log10(ns), log10(F(c, :)), 1) * [1; 0], 1:8));
            tc.verifyEqual(p(dma_curva(tc.R, ns)), 0.5, 'AbsTol', 0.1);
            tc.verifyEqual(p(dma_curva(tc.B, ns)), 1.5, 'AbsTol', 0.15);
            ms = [4 8 16 32 64];
            [Gv, Ga] = agregados_curva(tc.R, ms);
            q = @(G) mean(arrayfun(@(c) polyfit(log10(ms), log10(G(c, :)), 1) * [1; 0], 1:8));
            tc.verifyEqual(q(Gv), -1, 'AbsTol', 0.15);              % 2H - 2 con H = 0.5
            tc.verifyEqual(q(Ga), -0.5, 'AbsTol', 0.1);             % H - 1
        end

        function waveletYPsd(tc)
            E = wavelet_curva(tc.R, 'db4', 6);
            j = 1:5;
            p = @(E) mean(arrayfun(@(c) polyfit(j, log2(E(c, j)), 1) * [1; 0], 1:8));
            tc.verifyEqual(p(E), 0, 'AbsTol', 0.2);                  % ruido blanco: 2H - 1 = 0
            tc.verifyEqual(p(wavelet_curva(tc.B, 'db4', 6)), 2, 'AbsTol', 0.3);
            [P, f] = psd_curva(tc.R, 250, 256);
            b = f >= 2 & f <= 60;
            s = mean(arrayfun(@(c) polyfit(log10(f(b)), log10(P(c, b)), 1) * [1; 0], 1:8));
            tc.verifyEqual(s, 0, 'AbsTol', 0.2);
            [P, f] = psd_curva(tc.B, 250, 256);
            s = mean(arrayfun(@(c) polyfit(log10(f(b)), log10(P(c, b)), 1) * [1; 0], 1:8));
            tc.verifyEqual(s, -2, 'AbsTol', 0.3);
        end

        function boxcountRango(tc)
            ss = [64 128 256 512 1024];                          % escalas finas: la senoide es lisa
            pend = @(V) mean(arrayfun(@(c) polyfit(log10(ss), log10(boxcount_curva(V(c, :), ss)), 1) * [1; 0], 1:8));
            dR = pend(tc.R); dS = pend(tc.S);
            tc.verifyGreaterThan(dR, dS);
            tc.verifyGreaterThan(dR, 1.5); tc.verifyLessThan(dR, 2.05);
            tc.verifyLessThan(dS, 1.3);
        end

        function dimensionesDirectas(tc)
            DR = dimensiones_directas(tc.R); DS = dimensiones_directas(tc.S);
            tc.verifyTrue(all(mean(DR(:, 1:4)) > mean(DS(:, 1:4))));   % Katz, Petrosian, Sevcik
            f = (1:8).' / 256;
            tc.verifyEqual(DS(:, 5), 2 * sin(pi * f), 'RelTol', 0.02);  % movilidad de una senoide
            tc.verifyEqual(DS(:, 6), ones(8, 1), 'AbsTol', 0.05);       % complejidad de Hjorth = 1
        end

        function entropias(tc)
            x = tc.R(1:2, 1:600); s = tc.S(1:2, 1:600);
            SR = sampen_rejilla(x, 2, 0.2, 1);
            tc.verifyEqual(mean(SR), 2.2, 'AbsTol', 0.3);          % SampEn(2, 0.2) del ruido blanco
            tc.verifyLessThan(mean(sampen_rejilla(s, 2, 0.2, 1)), 0.5);
            S2 = sampen_rejilla(x, [1 2], [0.15 0.2], [1 2]);
            tc.verifySize(S2, [2 8]);
            tc.verifyEqual(S2(:, 4), SR, 'AbsTol', 1e-12);          % tau 1, m 2, r 0.2
            FR = fuzzyen_rejilla(x, 2, 0.2, 1); FS = fuzzyen_rejilla(s, 2, 0.2, 1);
            tc.verifyGreaterThan(mean(FR), mean(FS));
            P = permen_rejilla([x; s], [3 4], 1);
            tc.verifyGreaterThan(mean(P(1:2, 1)), 0.98);              % PE(3) del ruido ~ 1
            tc.verifyLessThan(mean(P(3:4, 1)), 0.6);
            L = lzc_rejilla([x; s], 1);
            tc.verifyGreaterThan(mean(L(1:2, 1)), 0.85);
            tc.verifyLessThan(mean(L(3:4, 1)), 0.3);
            M = mse_ci(x, [1 3], 2, 0.15);
            tc.verifyGreaterThan(M(:, 2), M(:, 1));
        end
    end
end
