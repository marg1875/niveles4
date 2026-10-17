function cls = envolver_seleccion(base, sel)
%ENVOLVER_SELECCION  Selección de características DENTRO de cada fold.
%   cls = envolver_seleccion(base, sel) envuelve el clasificador en
%   weka.classifiers.meta.AttributeSelectedClassifier: la selección se
%   recalcula con los datos de entrenamiento de cada fold, de modo que las
%   instancias de prueba nunca influyen en qué atributos se eligen.
%     'ninguna' : sin selección (devuelve base)
%     'ig<k>'   : ganancia de información + Ranker, mejores k atributos
%     'cfs'     : CFS (correlación atributo-clase vs redundancia) + BestFirst
if strcmp(sel, 'ninguna')
    cls = base; return
end
cls = weka.classifiers.meta.AttributeSelectedClassifier();
cls.setClassifier(base);
if startsWith(sel, 'ig')
    ev = weka.attributeSelection.InfoGainAttributeEval();
    bus = weka.attributeSelection.Ranker();
    bus.setNumToSelect(str2double(sel(3:end)));
elseif strcmp(sel, 'cfs')
    ev = weka.attributeSelection.CfsSubsetEval();
    bus = weka.attributeSelection.BestFirst();
else
    error('envolver_seleccion:sel', 'Selección desconocida: %s', sel);
end
cls.setEvaluator(ev);
cls.setSearch(bus);
end
