classdef PruebasAjuste < matlab.unittest.TestCase
%PRUEBASAJUSTE  Etapa 6: clasificadores parametrizados, selección, CV repetida y anidada.
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
                {'N00','N10','N40'}, y, {});
            tc.addTeardown(@() delete(tc.arff));
        end
    end

    methods (Test)
        function defectosIgualesAlArticulo(tc)
            % Sin parámetros, weka_clasificador_param = weka_clasificador
            % (valores de la CLI de Weka: 57.78 / 66.67 / 62.22 %)
            D = weka_cargar_arff(tc.arff);
            esp = {'kNN', 57.7778; 'SVM', 66.6667; 'MLP', 62.2222};
            for i = 1:3
                m = especificar_modelo(esp{i,1}, struct());
                M = weka_validar(D, weka_clasificador_param(m, tc.cfg), tc.cfg);
                tc.verifyEqual(M.acc, esp{i,2}, 'AbsTol', 1e-3, esp{i,1});
            end
        end

        function rejillasCompletas(tc)
            M = rejillas_hiperparametros(tc.cfg);
            R = tc.cfg.ajuste.rejilla; ns = numel(tc.cfg.ajuste.selecciones);
            nknn = numel(R.knn_K)*numel(R.knn_dist)*numel(R.knn_peso);
            nsvm = numel(R.svm_filtro)*(numel(R.svm_E)*numel(R.svm_C_poly) + numel(R.svm_gamma)*numel(R.svm_C_rbf));
            nmlp = numel(R.mlp_H)*numel(R.mlp_L)*numel(R.mlp_N)*numel(R.mlp_V);
            nrf  = numel(R.rf_arboles)*numel(R.rf_K)*numel(R.rf_prof);
            tc.verifyEqual(numel(M), ns*(nknn + nsvm + nmlp + nrf));
            tc.verifyEqual(numel(unique({M.id})), numel(M));        % ids únicos
        end

        function todosLosModelosCorren(tc)
            D = weka_cargar_arff(tc.arff);
            ms = [especificar_modelo('kNN', struct('K', 5, 'dist', 'euclid', 'peso', 'inv'), 'ig2'), ...
                  especificar_modelo('SVM', struct('kernel', 'rbf', 'gamma', 0.1, 'C', 10, 'filtro', 1), 'cfs'), ...
                  especificar_modelo('MLP', struct('H', '3', 'L', 0.1, 'N', 100, 'V', 20)), ...
                  especificar_modelo('RF', struct('arboles', 50, 'K', 2, 'prof', 5))];
            for m = ms
                M = weka_validar(D, weka_clasificador_param(m, tc.cfg), tc.cfg);
                tc.verifyGreaterThan(M.acc, 33, m.id);
            end
            v = especificar_modelo('Vote', struct('miembros', 'kNN+SVM'));
            v.miembros = ms(1:2);
            tc.verifyGreaterThan(weka_validar(D, weka_clasificador_param(v, tc.cfg), tc.cfg).acc, 33);
        end

        function seleccionDentroDelFold(tc)
            % El envoltorio entrena la selección con los datos que recibe:
            % con k = 2 el modelo resultante usa exactamente 2 atributos.
            D = weka_cargar_arff(tc.arff);
            cls = envolver_seleccion(weka.classifiers.lazy.IBk(), 'ig2');
            cls.buildClassifier(D);
            txt = char(cls.toString());
            k = regexp(txt, 'Selected attributes:[^:]*:\s*(\d+)', 'tokens', 'once');
            tc.verifyNotEmpty(k, txt);
            tc.verifyEqual(str2double(k{1}), 2);
        end

        function repetidaConUnaSemillaIgualASimple(tc)
            D = weka_cargar_arff(tc.arff);
            m = especificar_modelo('kNN', struct());
            Mm = weka_cv_repetida(D, m, tc.cfg, 1);
            tc.verifyEqual(Mm.acc, 57.7778, 'AbsTol', 1e-3);
            Mm = weka_cv_repetida(D, m, tc.cfg, 1:3);
            tc.verifyNumElements(Mm.acc_repeticiones, 3);
        end

        function anidadaUsaSoloEntrenamiento(tc)
            % Con un único candidato, la CV anidada se reduce a la CV externa
            % estándar de Weka (mismas particiones) -> mismo resultado que
            % crossValidateModel: la búsqueda interna no toca el fold de prueba.
            D = weka_cargar_arff(tc.arff);
            m = especificar_modelo('kNN', struct());
            [M, eleg] = weka_cv_anidada(D, m, tc.cfg);
            tc.verifyEqual(M.acc, 57.7778, 'AbsTol', 1e-3);
            tc.verifyEqual(M.n, 90);
            tc.verifyNumElements(eleg, 5);
            % Con varios candidatos devuelve uno elegido por fold
            c = [m, especificar_modelo('kNN', struct('K', 1)), especificar_modelo('SVM', struct())];
            [M2, eleg2] = weka_cv_anidada(D, c, tc.cfg);
            tc.verifyEqual(M2.n, 90);
            tc.verifyTrue(all(ismember(eleg2, {c.id})));
        end
    end
end
