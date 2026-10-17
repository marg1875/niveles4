function out = cv_anidada_conjuntos(ruta, cand, m, cfg)
%CV_ANIDADA_CONJUNTOS  CV anidada de la elección del conjunto de familias.
%   En cada fold externo se elige el conjunto (cand) con CV interna sobre el
%   entrenamiento y se evalúa en el fold de prueba con el modelo m.
%   out.acc (%), out.elegidos (conjunto elegido en cada fold).
weka_iniciar(cfg);
data = weka_cargar_arff(ruta);
k = cfg.weka.folds;
rnd = java.util.Random(cfg.weka.semilla);
D = weka.core.Instances(data); D.randomize(rnd); D.stratify(k);
atributos = @(cj) cellfun(@(f) strcat(f, '_', cfg.canales.nombres), strsplit(cj, '+'), 'UniformOutput', false);
aciertos = 0; total = 0; elegidos = cell(1, k);
for f = 0:k-1
    tr = D.trainCV(k, f, rnd); te = D.testCV(k, f);
    acc = zeros(1, numel(cand));
    for c = 1:numel(cand)
        a = atributos(cand{c}); a = [a{:}];
        acc(c) = weka_validar(weka_seleccionar_atributos(tr, a), weka_clasificador_param(m, cfg), cfg).acc;
    end
    [~, j] = max(acc);
    elegidos{f+1} = cand{j};
    a = atributos(cand{j}); a = [a{:}];
    trj = weka_seleccionar_atributos(tr, a); tej = weka_seleccionar_atributos(te, a);
    cls = weka_clasificador_param(m, cfg); cls.buildClassifier(trj);
    for i = 0:tej.numInstances()-1
        aciertos = aciertos + (cls.classifyInstance(tej.instance(i)) == tej.instance(i).classValue());
    end
    total = total + tej.numInstances();
end
out.acc = 100 * aciertos / total;
out.elegidos = elegidos;
end
