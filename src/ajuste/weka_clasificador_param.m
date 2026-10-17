function cls = weka_clasificador_param(m, cfg)
%WEKA_CLASIFICADOR_PARAM  Construye un clasificador de Weka a partir de una especificación.
%   cls = weka_clasificador_param(m, cfg), con m una especificación de
%   modelo (ver especificar_modelo):
%     m.clf    : 'kNN' | 'SVM' | 'MLP' | 'RF' | 'Vote'
%     m.params : struct con los hiperparámetros (los ausentes toman el
%                valor del artículo, ver weka_clasificador)
%     m.sel    : selección de características dentro de cada fold
%                ('ninguna' | 'ig8' | 'ig16' | 'ig24' | 'cfs')
%     m.miembros (solo Vote): arreglo de especificaciones
%
%   Hiperparámetros reconocidos:
%     kNN : K, dist ('euclid'|'manhattan'|'chebyshev'), peso ('ninguna'|'inv')
%     SVM : kernel ('poly'|'rbf'), E (poly), gamma (rbf), C, filtro (0 normaliza | 1 estandariza)
%     MLP : H (texto Weka: 'a', '7', '7,7'), L, M, N (épocas), V (% validación,
%           0 = sin parada temprana), D (decay: 0/1)
%     RF  : arboles, K (atributos por división, 0 = defecto), prof (0 = sin límite)
%     LOG : ridge (= 1/(2C) de sklearn; 0.5 por defecto), iter
p = m.params;
switch m.clf
    case 'kNN'
        cls = weka.classifiers.lazy.IBk();
        cls.setKNN(obt(p, 'K', 10));
        nn = weka.core.neighboursearch.LinearNNSearch();
        switch obt(p, 'dist', 'chebyshev')
            case 'euclid',    nn.setDistanceFunction(weka.core.EuclideanDistance());
            case 'manhattan', nn.setDistanceFunction(weka.core.ManhattanDistance());
            case 'chebyshev', nn.setDistanceFunction(weka.core.ChebyshevDistance());
        end
        cls.setNearestNeighbourSearchAlgorithm(nn);
        if strcmp(obt(p, 'peso', 'ninguna'), 'inv')
            cls.setDistanceWeighting(weka.core.SelectedTag( ...
                weka.classifiers.lazy.IBk.WEIGHT_INVERSE, weka.classifiers.lazy.IBk.TAGS_WEIGHTING));
        end
    case 'SVM'
        cls = weka.classifiers.functions.SMO();
        cls.setC(obt(p, 'C', 1));
        if strcmp(obt(p, 'kernel', 'poly'), 'rbf')
            k = weka.classifiers.functions.supportVector.RBFKernel();
            k.setGamma(obt(p, 'gamma', 0.01));
        else
            k = weka.classifiers.functions.supportVector.PolyKernel();
            k.setExponent(obt(p, 'E', cfg.weka.svm_exponente));
        end
        cls.setKernel(k);
        filtro = weka.classifiers.functions.SMO.FILTER_NORMALIZE;
        if obt(p, 'filtro', 0) == 1, filtro = weka.classifiers.functions.SMO.FILTER_STANDARDIZE; end
        cls.setFilterType(weka.core.SelectedTag(filtro, weka.classifiers.functions.SMO.TAGS_FILTER));
    case 'MLP'
        cls = weka.classifiers.functions.MultilayerPerceptron();
        cls.setHiddenLayers(obt(p, 'H', '7'));
        cls.setLearningRate(obt(p, 'L', 0.3));
        cls.setMomentum(obt(p, 'M', 0.2));
        cls.setTrainingTime(obt(p, 'N', 500));
        cls.setValidationSetSize(obt(p, 'V', 0));
        cls.setValidationThreshold(20);
        cls.setDecay(logical(obt(p, 'D', 0)));
    case 'RF'
        cls = weka.classifiers.trees.RandomForest();
        cls.setNumIterations(obt(p, 'arboles', 100));
        cls.setNumFeatures(obt(p, 'K', 0));
        cls.setMaxDepth(obt(p, 'prof', 0));
        cls.setBreakTiesRandomly(true);
    case 'LOG'
        % Regresión logística (ganadora en Niveles3: sklearn LogisticRegression
        % con StandardScaler, L2, C = 1). Weka minimiza -logL + ridge*||b||^2 y
        % sklearn 0.5||w||^2 + C*sum(pérdida) -> ridge = 1/(2C).
        lg = weka.classifiers.functions.Logistic();
        lg.setRidge(obt(p, 'ridge', 0.5));
        lg.setMaxIts(int32(obt(p, 'iter', 2000)));
        cls = weka.classifiers.meta.FilteredClassifier();
        cls.setFilter(weka.filters.unsupervised.attribute.Standardize());
        cls.setClassifier(lg);
    case 'Vote'
        cls = weka.classifiers.meta.Vote();
        arr = javaArray('weka.classifiers.Classifier', numel(m.miembros));
        for i = 1:numel(m.miembros)
            arr(i) = weka_clasificador_param(m.miembros(i), cfg);
        end
        cls.setClassifiers(arr);
        cls.setCombinationRule(weka.core.SelectedTag( ...
            weka.classifiers.meta.Vote.AVERAGE_RULE, weka.classifiers.meta.Vote.TAGS_RULES));
        return                       % los miembros ya llevan su propia selección
    otherwise
        error('weka_clasificador_param:clf', 'Clasificador desconocido: %s', m.clf);
end
cls = envolver_seleccion(cls, m.sel);
end

function v = obt(p, campo, defecto)
if isfield(p, campo), v = p.(campo); else, v = defecto; end
end
