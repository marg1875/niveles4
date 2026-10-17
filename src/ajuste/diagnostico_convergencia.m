function T = diagnostico_convergencia(ruta, cfg)
%DIAGNOSTICO_CONVERGENCIA  Curvas de convergencia/sensibilidad en una sesión.
%   T = diagnostico_convergencia(ruta, cfg) evalúa sobre el ARFF de la
%   sesión (configuración y conjunto de cfg.ajuste):
%   - MLP: precisión de ENTRENAMIENTO (modelo ajustado y evaluado sobre todos
%     los datos) y de CV 5-fold para N épocas x tasa de aprendizaje L, con y
%     sin parada temprana (V = 20 % de validación). Entrenamiento ~100 % con
%     CV estancada = sobreajuste; entrenamiento bajo que no sube = no converge.
%   - RF: precisión de CV vs número de árboles (debe estabilizarse).
%   - Conjuntos de familias de Hurst con los clasificadores del artículo.
%   Devuelve una tabla larga: prueba, clasificador, parametro, valor, extra,
%   acc_cv, acc_entrenamiento.
D = cargar_datos_ajuste(ruta, cfg);
filas = {};
for L = cfg.ajuste.diag_mlp_L
    for V = [0 20]
        for N = cfg.ajuste.diag_mlp_N
            m = especificar_modelo('MLP', struct('H', '7', 'L', L, 'N', N, 'V', V));
            cv = weka_validar(D, weka_clasificador_param(m, cfg), cfg);
            cls = weka_clasificador_param(m, cfg);
            cls.buildClassifier(D);
            ev = weka.classifiers.Evaluation(D);
            ev.evaluateModel(cls, D, nulo(D));
            filas(end+1,:) = {'mlp_epocas', 'MLP', sprintf('L=%g;V=%d', L, V), N, NaN, cv.acc, ev.pctCorrect()}; %#ok<AGROW>
        end
    end
end
for a = cfg.ajuste.diag_rf_arboles
    m = especificar_modelo('RF', struct('arboles', a));
    cv = weka_validar(D, weka_clasificador_param(m, cfg), cfg);
    filas(end+1,:) = {'rf_arboles', 'RF', 'arboles', a, NaN, cv.acc, NaN}; %#ok<AGROW>
end
conjs = cfg.ajuste.diag_conjuntos;
for i = 1:numel(conjs)
    Di = cargar_datos_ajuste(ruta, cfg, conjs{i});
    for c = {'kNN', 'SVM', 'MLP'}
        cv = weka_validar(Di, c{1}, cfg);
        filas(end+1,:) = {'conjunto', c{1}, conjs{i}, NaN, Di.numAttributes()-1, cv.acc, NaN}; %#ok<AGROW>
    end
end
T = cell2table(filas, 'VariableNames', {'prueba','clasificador','parametro','valor', ...
    'n_atributos','acc_cv','acc_entrenamiento'});
end

function extra = nulo(D)
salida = weka.classifiers.evaluation.output.prediction.Null();
salida.setBuffer(java.lang.StringBuffer());
salida.setHeader(D);
extra = javaArray('java.lang.Object', 1);
extra(1) = salida;
end
