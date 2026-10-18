classdef PruebasFractales < matlab.unittest.TestCase
%PRUEBASFRACTALES  Algoritmos fractales de Niveles3 (RS, Higuchi, DFA, VAR).
%   Valores teóricos: ruido blanco (H = 0.5, D = 2, alfa = 0.5) y movimiento
%   browniano = suma acumulada del ruido (D = 1.5, alfa = 1.5, H_var = 0.5).
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
        function higuchiValoresTeoricos(tc)
            rng(11); R = randn(20, 2048);
            tc.verifyEqual(mean(higuchi_dimension(R, tc.cfg)), 2, 'AbsTol', 0.1);
            tc.verifyEqual(mean(higuchi_dimension(cumsum(R, 2), tc.cfg)), 1.5, 'AbsTol', 0.1);
            s = sin(2*pi*(0:2047)/512);
            tc.verifyLessThan(higuchi_dimension(s, tc.cfg), 1.1);
        end

        function dfaValoresTeoricos(tc)
            rng(12); R = randn(20, 2048);
            tc.verifyEqual(mean(dfa_exponente(R, tc.cfg)), 0.5, 'AbsTol', 0.1);
            tc.verifyEqual(mean(dfa_exponente(cumsum(R, 2), tc.cfg)), 1.5, 'AbsTol', 0.15);
        end

        function variogramaBrowniano(tc)
            rng(13); B = cumsum(randn(20, 2048), 2);
            tc.verifyEqual(mean(variograma_pendiente(B, tc.cfg)), 0.5, 'AbsTol', 0.1);
        end

        function rsRuidoYTendencia(tc)
            rng(14); R = randn(20, 2048);
            h = hurst_rs_tendencia(R, tc.cfg);
            tc.verifyGreaterThan(mean(h), 0.45); tc.verifyLessThan(mean(h), 0.65);
            ht = hurst_rs_tendencia(R + 0.01 * (1:2048), tc.cfg);   % tendencia lineal
            tc.verifyEqual(ht, h, 'AbsTol', 1e-10);
        end

        function nanEnSenalCortaOConstante(tc)
            c = tc.cfg;
            tc.verifyTrue(isnan(hurst_rs_tendencia(randn(1, 30), c)));
            tc.verifyTrue(isnan(dfa_exponente(ones(1, 512), c)));
            tc.verifyTrue(isnan(higuchi_dimension(ones(1, 512), c)));
            tc.verifyTrue(isnan(variograma_pendiente(ones(1, 512), c)));
        end

        function vectorizadoIgualCanalACanal(tc)
            rng(15); V = cumsum(randn(4, 512), 2) + randn(4, 512);
            fun = {@hurst_rs_tendencia, @higuchi_dimension, @dfa_exponente, @variograma_pendiente};
            for i = 1:numel(fun)
                todo = fun{i}(V, tc.cfg);
                uno = arrayfun(@(c) fun{i}(V(c, :), tc.cfg), (1:4).');
                tc.verifyEqual(todo, uno, 'AbsTol', 1e-12);
            end
        end

        function candidatosDesdeCurvasCoinciden(tc)
            rng(17); V = cumsum(randn(16, 512), 2) + 3 * randn(16, 512);
            c = tc.cfg;
            baseH = struct('nombre', 'Re16a2', 'modo', 'ecuacion13', 'm_max', 16, 'm_min', 2);
            [K, lay] = candidatos_hig_dfa(c, baseH);
            [X, n] = caracteristicas_desde_curvas(curvas_ventana(V, c, baseH, lay), K, lay, c);
            tc.verifyEqual(numel(n), numel(K) * 16);
            % HIG1a30 (curvas) = Higuchi por defecto de Niveles3 (HIGn3)
            tc.verifyEqual(X(startsWith(n, 'HIG1a30_')), X(startsWith(n, 'HIGn3_')), 'AbsTol', 1e-10);
            tc.verifyEqual(X(startsWith(n, 'HIG4a20_')).', higuchi_dimension(V, c, 4, 20), 'AbsTol', 1e-10);
            % DFA1 en un rango = dfa_exponente con esas escalas
            ns = c.calib_hd.dfa_ns(7:14);
            cc = c; cc.fractales.dfa_min_escalas = 4;
            tc.verifyEqual(X(startsWith(n, sprintf('DFAl%da%d_', ns(1), ns(end)))).', ...
                dfa_exponente(V, cc, ns, 1), 'AbsTol', 1e-10);
            tc.verifyEqual(X(startsWith(n, 'Re16a2_')), caracterizadores_limites(V, baseH, c), 'AbsTol', 1e-12);
        end

        function dfa2BrownianoYTendenciaCuadratica(tc)
            rng(18); B = cumsum(randn(10, 2048), 2);
            c = tc.cfg; ns = c.calib_hd.dfa_ns(5:end-1);
            a2 = dfa_exponente(B, c, ns, 2);
            tc.verifyEqual(mean(a2), 1.5, 'AbsTol', 0.15);
            t = linspace(-1, 1, 2048);
            R = randn(10, 2048);
            a = dfa_exponente(R + 5 * t.^2, c, ns, 2);      % DFA2 elimina tendencias cuadráticas del perfil
            tc.verifyEqual(mean(a), 0.5, 'AbsTol', 0.1);
        end

        function dfa3YTraslape(tc)
            rng(19); R = randn(10, 2048); B = cumsum(R, 2);
            c = tc.cfg; ns = [8 10 12 16 20 24 32 40 48 64 80 96 128];
            tc.verifyEqual(mean(dfa_exponente(R, c, ns, 3)), 0.5, 'AbsTol', 0.1);
            tc.verifyEqual(mean(dfa_exponente(B, c, ns, 3)), 1.5, 'AbsTol', 0.15);
            F1 = dfa_curva(R, ns, 1, 1); F2 = dfa_curva(R, ns, 1, 0.5);
            tc.verifyEqual(mean(F2), mean(F1), 'RelTol', 0.05);   % traslape: misma F (en promedio), menos varianza
            % rejilla por longitud de ventana: límites y nombres
            cc = c; cc.calib_hd.dfa_ordenes = 1:3; cc.calib_hd.dfa_pasos = [1 0.5]; cc.calib_hd.prefijo = 'W250';
            [K, lay] = candidatos_hig_dfa(cc, struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {}), 250);
            tc.verifyLessThanOrEqual(lay.kmax, 62);
            tc.verifyTrue(any(strcmp({K.nombre}, 'W250DFAto8a64')));
            tc.verifyFalse(any(strcmp({K.nombre}, 'W250DFAt8a128')));   % 128 > 250/3 sin traslape
            V = randn(16, 250);
            [X, n] = caracteristicas_desde_curvas(curvas_ventana(V, cc, struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {}), lay), K, lay, cc);
            tc.verifyEqual(numel(n), 16 * numel(K));
            tc.verifyEqual(X(startsWith(n, 'W250DFAto8a64_')).', dfa_exponente(V, struct('fractales', struct('dfa_min_escalas', 4)), ...
                lay.ns(lay.ns >= 8 & lay.ns <= 64), 3, 0.5), 'AbsTol', 1e-10);
        end

        function ventanaNombresYOrden(tc)
            rng(16); V = randn(16, 512);
            base = struct('nombre', 'Re16a2', 'modo', 'ecuacion13', 'm_max', 16, 'm_min', 2);
            [f, n] = fractales_ventana(V, tc.cfg, base);
            tc.verifySize(f, [1 16 * 9]);
            tc.verifyEqual(n{1}, 'Re16a2_Fp1');
            tc.verifyEqual(n{17}, 'RS_Fp1');
            i = find(strcmp(n, 'HIG_Cz'));
            tc.verifyEqual(f(i), higuchi_dimension(V(6, :), tc.cfg), 'AbsTol', 1e-12);
            ip = find(strcmp(n, 'PROM_Cz'));
            fam = cellfun(@(a) f(strcmp(n, [a '_Cz'])), {'RS','HIG','DFA','VAR'});
            tc.verifyEqual(f(ip), mean(fam), 'AbsTol', 1e-12);
            tc.verifyTrue(all(isfinite(f)));
        end
    end
end
