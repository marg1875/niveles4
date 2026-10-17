classdef PruebasHurst < matlab.unittest.TestCase
%PRUEBASHURST  Fidelidad y validez de HO / HRS / HV.
%   Oráculo: pruebas/referencia/rangors.m y variograma.m (copias sin
%   modificar del código del grupo, 1PM10/Analisis matlab).

    properties
        cfg
    end

    methods (TestClassSetup)
        function rutas(tc)
            raiz = fileparts(fileparts(mfilename('fullpath')));
            addpath(raiz, genpath(fullfile(raiz, 'src')), fullfile(raiz, 'pruebas', 'referencia'));
            tc.cfg = config_niveles();
            w = warning('off', 'MATLAB:rankDeficientMatrix');
            tc.addTeardown(@() warning(w));
        end
    end

    methods (Test)
        function nivelesIgualesARangors(tc)
            % Los tamaños de bloque de niveles_rangors = tau de rangors.m
            rng(1);
            for N = [25 50 75 100 125 250 375 500 512 700]
                v = randn(1, N);
                for p = [2 4 8 9 16 32 64 128]
                    if N/p < 2, continue, end
                    tau = rangors(v, p);
                    tc.verifyEqual(niveles_rangors(N, p, 1), tau, sprintf('N=%d p=%d', N, p));
                end
            end
        end

        function rsIgualAOriginalModoRangors(tc)
            % En modo 'rangors_original' se reproduce exactamente rangors.m
            % (que usa '/' matricial, no la media de R./S).
            rng(2);
            for N = [100 512 700]
                v = randn(1, N);
                for p = [4 9 64]
                    [tau, RS] = rangors(v, p);
                    ok = tau >= 2;
                    tc.verifyEqual(rs_diadico(v, tau(ok), 'rangors_original'), RS(ok), 'RelTol', 1e-12);
                end
            end
        end

        function rsIgualAEcuacion13(tc)
            % Modo por defecto = media de R./S (Ec. 13 del artículo)
            rng(3); v = randn(1, 512);
            for m = [512 256 64 8]
                B = reshape(v, m, []);
                Y = cumsum(B - mean(B, 1), 1);
                esperado = mean((max(Y) - min(Y)) ./ std(B, 0, 1));
                tc.verifyEqual(rs_diadico(v, m), esperado, 'RelTol', 1e-12);
            end
        end

        function semivariogramaIgualABucle(tc)
            rng(4); v = randn(3, 300);
            G = semivariograma(v);
            for c = 1:3
                for h = [1 2 10 149 150]
                    esp = sum((v(c,1+h:end) - v(c,1:end-h)).^2) / (2*(300-h));
                    tc.verifyEqual(G(c,h), esp, 'AbsTol', 1e-10);
                end
            end
            % Coincide con variograma.m del grupo en los lags 2..N/2-1 (el
            % original no calcula el lag 1 y deja vario(1) = 0)
            vr = variograma(v(1,:), 300);
            tc.verifyEqual(G(1, 2:numel(vr)), vr(2:end), 'AbsTol', 1e-10);
        end

        function ruidoBlancoHOcercaDeMedio(tc)
            % R/S de ruido blanco con N = 512: H ~ 0.55-0.65 (sesgo de
            % Anis-Lloyd en muestras cortas)
            rng(5);
            H = hurst_ho(randn(200, 512), tc.cfg);
            tc.verifyGreaterThan(mean(H), 0.5);
            tc.verifyLessThan(mean(H), 0.7);
        end

        function monotonoEnH(tc)
            % fGn sintético (incrementos de wfbm): H estimado crece con H real
            if isempty(which('wfbm')), tc.assumeFail('Sin Wavelet Toolbox'); end
            rng(6);
            Hs = [0.3 0.5 0.7 0.85];
            est = zeros(size(Hs));
            for i = 1:numel(Hs)
                x = zeros(40, 512);
                for r = 1:40, b = wfbm(Hs(i), 513); x(r,:) = diff(b(:).'); end
                est(i) = mean(hurst_ho(x, tc.cfg));
            end
            tc.verifyTrue(all(diff(est) > 0), mat2str(est, 3));
        end

        function validezVentanasCortas(tc)
            % 25 muestras: HO necesita >= 3 escalas con bloque >= 8 -> NaN
            [f, ne] = extraer_hurst_ventana(randn(16, 25), tc.cfg);
            nombres = nombres_caracteristicas(tc.cfg);
            tc.verifyTrue(all(isnan(f(startsWith(nombres, 'HO_')))));
            tc.verifyEqual(ne(1), 2);
            % 512 muestras: todo estimable
            f = extraer_hurst_ventana(randn(16, 512), tc.cfg);
            tc.verifyFalse(any(isnan(f)));
            tc.verifyEqual(numel(f), numel(nombres));
        end

        function ventanaIgualAFuncionesSueltas(tc)
            rng(7); V = randn(4, 512);
            c = tc.cfg; c.canales.nombres = {'A','B','C','D'};
            f = extraer_hurst_ventana(V, c);
            n = nombres_caracteristicas(c);
            tc.verifyEqual(f(startsWith(n, 'HO_')).', hurst_ho(V, c), 'AbsTol', 1e-12);
            tc.verifyEqual(f(startsWith(n, 'HRS64_')).', hurst_hrs(V, 64, c), 'AbsTol', 1e-12);
            tc.verifyEqual(f(startsWith(n, 'HV_')).', hurst_hv(V, c), 'AbsTol', 1e-12);
        end
    end
end
