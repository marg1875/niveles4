classdef PruebasEventos < matlab.unittest.TestCase
%PRUEBASEVENTOS  Detección de ERD, eventos FIS y ventanas de reposo (etapa 8).
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
        function erdRecuperaEventosSinteticos(tc)
            % Ritmo de 12 Hz en los canales centro-parietales que se apaga
            % (ERD) durante 2 s en instantes conocidos.
            fs = tc.cfg.fs; N = 60*fs; t = (0:N-1)/fs;
            rng(1);
            X = 0.3*randn(16, N);
            amp = ones(1, N);
            centros = [12 24 36 48];
            for c = centros, amp(abs(t - c) < 1) = 0.1; end
            idx = ismember(tc.cfg.canales.nombres, tc.cfg.eventos.canales_erd);
            X(idx, :) = X(idx, :) + amp .* sin(2*pi*12*t);
            E = detectar_erd(X, tc.cfg);
            tc.verifyEqual(height(E), numel(centros));
            tc.verifyLessThanOrEqual(max(abs(E.minimo.'/fs - centros)), 0.5);
        end

        function erdIgnoraReposoInicial(tc)
            fs = tc.cfg.fs; N = 30*fs; t = (0:N-1)/fs;
            X = zeros(16, N);
            idx = ismember(tc.cfg.canales.nombres, tc.cfg.eventos.canales_erd);
            amp = ones(1, N); amp(t > 1 & t < 3) = 0.05;           % caída dentro de los 5 s
            X(idx, :) = repmat(amp .* sin(2*pi*12*t), nnz(idx), 1);
            E = detectar_erd(X, tc.cfg);
            tc.verifyTrue(all(E.inicio > tc.cfg.eventos.reposo_inicial_s * fs));
        end

        function eventosFisPorBloques(tc)
            st = repmat([0 0 0 0 40 40 40 40], 1, 10);
            E = eventos_fis(39000, st);
            tc.verifyEqual(height(E), 10);
            bordes = round((0:80) * 39000 / 80);
            tc.verifyEqual(E.inicio(1), bordes(5) + 1);
            tc.verifyEqual(E.fin(1), bordes(9));
        end

        function ventanasDefinicionB(tc)
            fs = tc.cfg.fs;
            E = table([3000; 8000], [4000; 9000], [3500; 8500], [NaN; NaN], ...
                'VariableNames', {'inicio','fin','minimo','profundidad'});
            W = ventanas_por_definicion(E, 20000, 'B', 'KMI', tc.cfg);
            R = W(strcmp(W.tipo, 'reposo'), :);
            tc.verifyEqual(height(R), 2);
            tc.verifyEqual(R.inicio, E.fin + round(0.5*fs));
            tc.verifyEqual(R.inicio + R.L - 1 - E.fin, repmat(round(0.5*fs) + tc.cfg.eventos.L_B - 1, 2, 1));
            tc.verifyLessThanOrEqual((R.L(1) + round(0.5*fs)) / fs, 0.81);   % termina hacia 0.8 s
            T = W(strcmp(W.tipo, 'tarea'), :);
            tc.verifyEqual(T.inicio + floor(T.L/2), E.minimo);               % centrada en el ERD
        end

        function ventanasDefinicionA(tc)
            fs = tc.cfg.fs;
            E = eventos_fis(39000, repmat([0 0 0 0 10 10 10 10], 1, 10));
            W = ventanas_por_definicion(E, 39000, 'A', 'FIS', tc.cfg);
            R = W(strcmp(W.tipo, 'reposo'), :);
            tc.verifyEqual(height(R), 4);
            tc.verifyLessThanOrEqual(max(R.inicio + R.L - 1), tc.cfg.eventos.reposo_inicial_s * fs);
            tc.verifyEqual(nnz(strcmp(W.tipo, 'tarea')), 30);
        end

        function caracterizadoresPorLimites(tc)
            % Con L = 512 los límites reproducen los rangos calibrados
            rng(4); V = randn(16, 512);
            fams = tc.cfg.eventos.familias;
            f = caracterizadores_limites(V, fams, tc.cfg);
            ms = floor(512 ./ 2.^(0:9)); ms = ms(ms >= 2);
            RSo = rs_diadico(V, ms, 'rangors_original');
            sel = ms <= 64;
            tc.verifyEqual(f(1:16).', pendiente_hurst(ms(sel), RSo(:, sel), 3), 'AbsTol', 1e-12);
            tc.verifyEqual(f(33:48).', hurst_ho(V, tc.cfg), 'AbsTol', 1e-12);   % HO
            % Con L = 75 sigue siendo estimable
            f75 = caracterizadores_limites(randn(16, 75), fams, tc.cfg);
            tc.verifyFalse(any(isnan(f75)));
        end
    end
end
