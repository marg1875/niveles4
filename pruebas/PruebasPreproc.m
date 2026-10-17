classdef PruebasPreproc < matlab.unittest.TestCase
%PRUEBASPREPROC  Filtros, blanqueo e ICA.
    properties
        cfg
    end

    methods (TestClassSetup)
        function rutas(tc)
            raiz = fileparts(fileparts(mfilename('fullpath')));
            addpath(raiz, genpath(fullfile(raiz, 'src')));
            tc.cfg = config_niveles();
        end
    end

    methods (Test)
        function notchAtenua60Hz(tc)
            F = disenar_filtros(tc.cfg);
            n = 10000; t = (0:n-1)/tc.cfg.fs;
            x = 5*sin(2*pi*60*t) + sin(2*pi*10*t);
            y = filtrar_eeg(x, F);
            X = abs(fft(x)); Y = abs(fft(y));
            i60 = round(60/tc.cfg.fs*n) + 1; i10 = round(10/tc.cfg.fs*n) + 1;
            tc.verifyGreaterThan(20*log10(X(i60)/Y(i60)), 30);
            tc.verifyEqual(Y(i10)/X(i10), 1, 'AbsTol', 0.02);   % banda de paso
        end

        function filtroFIRopcional(tc)
            c = tc.cfg; c.filtro.tipo = 'fir';
            F = disenar_filtros(c);
            y = filtrar_eeg(randn(2, 3000), F);
            tc.verifySize(y, [2 3000]);
            tc.verifyTrue(all(isfinite(y(:))));
        end

        function blanqueoDaIdentidad(tc)
            rng(1); X = randn(6) * randn(6, 5000);
            [Z, Q] = blanquear(X, tc.cfg);
            tc.verifyEqual(Z*Z.'/size(Z,2), eye(6), 'AbsTol', 1e-10);
            tc.verifySize(Q, [6 6]);
        end

        function blanqueoDescartaRangoDeficiente(tc)
            rng(2); X = randn(5, 4000); X(6,:) = X(1,:) + X(2,:);
            [~, Q] = blanquear(X, tc.cfg);
            tc.verifySize(Q, [5 6]);
        end

        function icaSeparaMezcla(tc)
            rng(7); n = 20000;
            S = [sign(randn(1,n)).*rand(1,n); rand(1,n)-0.5; trnd_simple(3,n); randn(1,n).^3];
            S = S ./ std(S, 0, 2);
            A = randn(4);
            [Z, Q] = blanquear(A*S, tc.cfg);
            [W, ~, info] = ica_gradiente_natural(Z, tc.cfg);
            P = abs(W*Q*A); k = 4;
            amari = (sum(sum(P./max(P,[],2),2)-1) + sum(sum(P./max(P,[],1),1)-1)) / (2*k*(k-1));
            tc.verifyLessThan(amari, 0.1);
            tc.verifyTrue(info.convergio);
            W2 = ica_gradiente_natural(Z, tc.cfg);
            tc.verifyEqual(W2, W);                                % determinista
        end

        function carEliminaModoComun(tc)
            rng(5); comun = 50 * randn(1, 5000);
            X = randn(16, 5000) + comun;                       % señal común dominante
            R = corrcoef(X.'); tc.verifyGreaterThan(min(R(:)), 0.99);
            Y = referencia_promedio(X);
            tc.verifyEqual(mean(Y, 1), zeros(1, 5000), 'AbsTol', 1e-10);
            R = corrcoef(Y.'); R(1:17:end) = 0;
            tc.verifyLessThan(max(abs(R(:))), 0.2);
        end

        function artefactosMarcadosYVentanas(tc)
            rng(6); X = randn(16, 10000);
            X(3, 5000:5010) = 40;                              % artefacto breve
            c = tc.cfg; c.artefactos.activar = true;          % desactivado por defecto (sensibilidad)
            malas = marcar_artefactos(X, c);
            tc.verifyEqual(size(malas, 1), 1);
            m = round(tc.cfg.artefactos.margen_s * tc.cfg.fs);
            tc.verifyLessThanOrEqual(malas(1, 1), 5000 - m + 1);
            tc.verifyGreaterThanOrEqual(malas(1, 2), 5010 + m - 1);
            tc.verifyEqual(ventana_valida(malas, [1 4800 5005 9000], 256), [true false false true]);
        end

        function puentesDetectadosYCorregidos(tc)
            rng(7); X = randn(16, 5000);
            X(4, :) = X(3, :) + 0.01 * randn(1, 5000);         % F3 y Fz puenteados
            c = tc.cfg; c.puentes.metodo = 'eliminar';
            [Y, p] = corregir_puentes(X, c);
            tc.verifyEqual(p, [3 4]);
            tc.verifyTrue(all(isnan(Y([3 4], :)), 'all'));
            c.puentes.metodo = 'interpolar';
            Y = corregir_puentes(X, c);
            tc.verifyTrue(all(isfinite(Y(:))));
            R = corrcoef(Y(3, :), Y(4, :));
            tc.verifyLessThan(abs(R(1, 2)), 0.999);
        end

        function rechazoMinimoDeCanales(tc)
            rng(8); X = randn(16, 10000);
            X(3, 2000:2010) = 40;                              % un solo canal
            c = tc.cfg; c.artefactos.activar = true; c.artefactos.min_canales = 2;
            tc.verifyEmpty(marcar_artefactos(X, c));
            X(5, 2000:2010) = 40;                              % ahora dos canales
            tc.verifyNotEmpty(marcar_artefactos(X, c));
        end

        function asrReparaRafaga(tc)
            rng(9); fs = tc.cfg.fs; n = 60 * fs;
            A = randn(16, 6);
            X = A * randn(6, n) + 0.3 * randn(16, n);          % EEG simulado (16 canales correlacionados)
            r = 30*fs+1 : 31*fs;                               % ráfaga de 1 s en 3 canales
            Xa = X; Xa(1:3, r) = Xa(1:3, r) + 40 * randn(3, numel(r));
            c = tc.cfg;
            E = asr_calibrar(Xa, c);
            [Y, frac] = asr_procesar(Xa, E, c);
            ruido = Xa(1:3, r) - X(1:3, r);
            resto = Y(1:3, r) - X(1:3, r);
            tc.verifyLessThan(norm(resto(:)) / norm(ruido(:)), 0.2);   % > 80 % de la ráfaga removida
            fuera = [1:25*fs, 36*fs:n];
            R = corrcoef(Y(:, fuera), X(:, fuera));
            tc.verifyGreaterThan(R(1, 2), 0.99);
            tc.verifyLessThan(frac, 0.15);
        end

        function asrNoTocaSenalLimpia(tc)
            rng(10); n = 60 * tc.cfg.fs;
            X = randn(16, 6) * randn(6, n) + 0.3 * randn(16, n);
            E = asr_calibrar(X, tc.cfg);
            Y = asr_procesar(X, E, tc.cfg);
            tc.verifyLessThan(max(abs(Y(:) - X(:))) / max(abs(X(:))), 0.05);
        end

        function curtosisDetectaArtefacto(tc)
            rng(3); Y = randn(8, 20000);
            Y(5, 1:50:end) = 40;                                  % "parpadeos"
            idx = limpiar_artefactos_ica(Y, tc.cfg);
            tc.verifyEqual(idx, 5);
        end
    end
end

function x = trnd_simple(nu, n)
% t de Student sin Statistics Toolbox: normal / sqrt(chi2/nu)
x = randn(1, n) ./ sqrt(sum(randn(nu, n).^2, 1) / nu);
end
