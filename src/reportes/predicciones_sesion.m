function P = predicciones_sesion(ruta, conjunto, m, cfg, semillas)
%PREDICCIONES_SESION  CV de 5 particiones de un ARFF guardando cada predicción.
%   P = predicciones_sesion(ruta, conjunto, m, cfg, semillas): mismas
%   particiones que crossValidateModel de Weka (randomize + stratify +
%   trainCV/testCV) con cada semilla; cada predicción se cruza con su fila
%   del *_meta.csv (archivo, clase, bloque). La usan analisis_confundido_reposo
%   y corrida_controles.
weka_iniciar(cfg);
data = weka_cargar_arff(ruta);
meta = readtable(strrep(ruta, '.arff', '_meta.csv'), 'TextType', 'char');
fam = strsplit(conjunto, '+');
attrs = {};
for f = 1:numel(fam), attrs = [attrs, strcat(fam{f}, '_', cfg.canales.nombres)]; end %#ok<AGROW>
D0 = weka_seleccionar_atributos(data, attrs);
% Atributo ID (primera posición) que el clasificador no verá
addid = weka.filters.unsupervised.attribute.AddID();
addid.setInputFormat(D0);
D0 = weka.filters.Filter.useFilter(D0, addid);
D0.setClassIndex(D0.numAttributes() - 1);
k = cfg.weka.folds;
filas = {};
for s = semillas
    rnd = java.util.Random(s);
    D = weka.core.Instances(D0);
    D.randomize(rnd);
    D.stratify(k);
    for f = 0:k-1
        tr = D.trainCV(k, f, rnd);
        te = D.testCV(k, f);
        fc = weka.classifiers.meta.FilteredClassifier();
        rm = weka.filters.unsupervised.attribute.Remove();
        rm.setAttributeIndices('first');
        fc.setFilter(rm);
        fc.setClassifier(weka_clasificador_param(m, cfg));
        fc.buildClassifier(tr);
        for i = 0:te.numInstances()-1
            inst = te.instance(i);
            fila = round(inst.value(0));                 % ID 1-based = fila del meta
            pred = char(te.classAttribute().value(fc.classifyInstance(inst)));
            filas(end+1, :) = {meta.sesion{fila}, meta.archivo{fila}, meta.nivel_archivo(fila), ...
                meta.clase(fila), meta.bloque(fila), meta.ancla(fila), s, pred}; %#ok<AGROW>
        end
    end
end
P = cell2table(filas, 'VariableNames', {'sesion','archivo','nivel_archivo','clase','bloque','ancla','semilla','prediccion'});
end
