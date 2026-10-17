classdef PruebasWeka < matlab.unittest.TestCase
%PRUEBASWEKA  ARFF y validación cruzada en Weka (Java API desde MATLAB).
    properties
        cfg
        arff
    end

    methods (TestClassSetup)
        function preparar(tc)
            raiz = fileparts(fileparts(mfilename('fullpath')));
            addpath(raiz, genpath(fullfile(raiz, 'src')));
            tc.cfg = config_niveles();
            tc.assumeNotEmpty(tc.cfg.rutas.weka_jar, 'weka.jar no encontrado');
            weka_iniciar(tc.cfg);
            rng(2); n = 90;
            X = randn(n, 4) + repelem([0; 0.8; 1.6], 30, 1);
            X(3, 2) = NaN;
            y = repelem({'N00'; 'N10'; 'N40'}, 30, 1);
            tc.arff = [tempname '.arff'];
            exportar_arff(tc.arff, 'prueba', {'HO_F3','HV_F3','HO_F4','HV_F4'}, X, ...
                {'N00','N10','N40'}, y, {'comentario de prueba'});
            tc.addTeardown(@() delete(tc.arff));
        end
    end

    methods (Test)
        function idaYVueltaARFF(tc)
            D = weka_cargar_arff(tc.arff);
            tc.verifyEqual(D.numInstances(), 90);
            tc.verifyEqual(D.numAttributes(), 5);
            tc.verifyEqual(D.attributeStats(1).missingCount, 1);   % el NaN -> '?'
            tc.verifyEqual(char(D.classAttribute().value(2)), 'N40');
        end

        function seleccionDeAtributos(tc)
            D = weka_cargar_arff(tc.arff);
            [S, ne] = weka_seleccionar_atributos(D, {'HO_F3', 'HO_F4'});
            tc.verifyEqual(S.numAttributes(), 3);
            tc.verifyEqual(ne, 2);
            tc.verifyEqual(char(S.attribute(0).name()), 'HO_F3');
        end

        function cvIgualAWekaCLI(tc)
            % Valores obtenidos con la CLI de Weka 3.8.6 sobre el mismo ARFF:
            %   java -cp weka.jar weka.classifiers.lazy.IBk -K 10 -A "...Chebyshev..." -x 5 -s 1
            %   -> 57.7778 %, kappa 0.3667, MAE 0.3345
            %   SMO -> 66.6667 % ; MLP -L 0.3 -M 0.2 -N 500 -H 7 -> 62.2222 %
            D = weka_cargar_arff(tc.arff);
            M = weka_validar(D, 'kNN', tc.cfg);
            tc.verifyEqual(M.acc, 57.7778, 'AbsTol', 1e-3);
            tc.verifyEqual(M.kappa, 0.3667, 'AbsTol', 1e-3);
            tc.verifyEqual(M.mae, 0.3345, 'AbsTol', 1e-3);
            tc.verifyEqual(weka_validar(D, 'SVM', tc.cfg).acc, 66.6667, 'AbsTol', 1e-3);
            tc.verifyEqual(weka_validar(D, 'MLP', tc.cfg).acc, 62.2222, 'AbsTol', 1e-3);
            tc.verifySize(M.confusion, [3 3]);
            tc.verifyEqual(sum(M.confusion(:)), 90);
        end

        function cvAgrupadaNoMezclaGrupos(tc)
            D = weka_cargar_arff(tc.arff);
            c = tc.cfg; c.weka.cv_modo = 'agrupada';
            M = weka_validar(D, 'kNN', c, repelem((1:30).', 3));
            tc.verifyEqual(M.n, 90);
        end

        function todosLosClasificadoresCorren(tc)
            D = weka_cargar_arff(tc.arff);
            for c = {'BN', 'NB', 'RF', 'RT'}
                M = weka_validar(D, c{1}, tc.cfg);
                tc.verifyGreaterThan(M.acc, 33, c{1});
            end
        end
    end
end
