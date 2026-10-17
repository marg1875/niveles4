classdef PruebasSegmentacion < matlab.unittest.TestCase
%PRUEBASSEGMENTACION  Bloques, anclas y ventanas.
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
        function bloquesYAnclas(tc)
            st = repmat([0 0 0 0 40 40 40 40], 1, 10);
            N = 39100;
            A = segmentar_por_estimulos(N, st, tc.cfg);
            tc.verifyEqual(numel(unique(A.bloque)), 20);
            tc.verifyEqual(height(A), 60);                        % 3 anclas x 20 bloques
            tc.verifyEqual(nnz(A.tarea), 30);
            % Bloques impares = reposo, pares = tarea
            tc.verifyEqual(A.tarea, mod(A.bloque, 2) == 0);
            % Ventanas de 512 dentro del bloque y respetando los márgenes
            bordes = round((0:80) * N / 80);
            m0 = round(tc.cfg.segmentacion.margen_inicio_s * tc.cfg.fs);
            for i = 1:height(A)
                b = A.bloque(i);
                ini = bordes(4*(b-1) + 1) + 1 + m0;
                fin = bordes(4*b + 1);
                tc.verifyGreaterThanOrEqual(A.centro(i) - 256, ini);
                tc.verifyLessThanOrEqual(A.centro(i) + 255, fin);
            end
        end

        function mismasAnclasEntreDuraciones(tc)
            st = repmat([0 0 0 0 10 10 10 10], 1, 10);
            X = repmat(1:39000, 2, 1);
            A = segmentar_por_estimulos(39000, st, tc.cfg);
            for L = [512 500 250 25]
                V = extraer_ventanas(X, A.centro, L);
                tc.verifySize(V, [2 L height(A)]);
                % la ventana está centrada en el ancla
                tc.verifyLessThanOrEqual(abs(mean(V(1,:,1)) - A.centro(1)), 1);
            end
        end
    end
end
