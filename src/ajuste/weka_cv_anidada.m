function [M, elegidos] = weka_cv_anidada(data, candidatos, cfg)
%WEKA_CV_ANIDADA  CV anidada: búsqueda interna de hiperparámetros por fold.
%   [M, elegidos] = weka_cv_anidada(data, candidatos, cfg)
%   Bucle externo: k folds estratificados (mismo procedimiento que
%   crossValidateModel de Weka: randomize(Random(semilla)) + stratify(k) +
%   trainCV/testCV). Para cada fold externo, SOLO con su conjunto de
%   entrenamiento, se evalúa cada especificación de 'candidatos' con una CV
%   interna de k folds y se elige la de mayor precisión (empate -> la de
%   menor índice). El modelo elegido se entrena con todo el entrenamiento
%   externo y se evalúa en el fold externo de prueba, que nunca intervino en
%   la elección. Las métricas se acumulan en un único objeto Evaluation.
%   elegidos{f} = id del modelo elegido en el fold externo f.
k = cfg.weka.folds;
rnd = java.util.Random(cfg.weka.semilla);
D = weka.core.Instances(data);
D.randomize(rnd);
D.stratify(k);
ev = weka.classifiers.Evaluation(D);
salida = weka.classifiers.evaluation.output.prediction.Null();
salida.setBuffer(java.lang.StringBuffer());
salida.setHeader(D);
extra = javaArray('java.lang.Object', 1);
extra(1) = salida;
elegidos = cell(1, k);
for f = 0:k-1
    tr = D.trainCV(k, f, rnd);
    te = D.testCV(k, f);
    acc = -inf(1, numel(candidatos));
    for j = 1:numel(candidatos)
        c = cfg; c.weka.cv_modo = 'estratificada';
        Mi = weka_validar(tr, weka_clasificador_param(candidatos(j), c), c);   % CV interna
        acc(j) = Mi.acc;
    end
    [~, jmej] = max(acc);
    elegidos{f+1} = candidatos(jmej).id;
    cls = weka_clasificador_param(candidatos(jmej), cfg);
    cls.buildClassifier(tr);
    ev.setPriors(tr);
    ev.evaluateModel(cls, te, extra);
end
nc = D.numClasses();
M.n = ev.numInstances();
M.acc = ev.pctCorrect();
M.kappa = ev.kappa();
M.mae = ev.meanAbsoluteError();
M.wf = ev.weightedFMeasure();
M.tpr = arrayfun(@(i) ev.truePositiveRate(i), 0:nc-1);
M.confusion = double(ev.confusionMatrix());
end
