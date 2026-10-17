%% ANALISIS_CONFUNDIDO_REPOSO  ¿El clasificador reconoce el reposo o la grabación?
%
%   La clase N00 (reposo) se toma de los bloques de reposo de los DOS
%   archivos de la sesión (el de 10 % y el de 40 %), que son grabaciones
%   distintas. Hipótesis de confundido: si el modelo aprende "de qué
%   archivo viene la ventana", el reposo del archivo de 10 % se predice
%   como N10 y el del archivo de 40 % como N40.
%
%   Se repite la CV 5-fold de Weka (mismas particiones que
%   crossValidateModel: randomize(Random(s)) + stratify(5) + trainCV/testCV)
%   con un atributo de identificación que el clasificador NO usa
%   (FilteredClassifier + Remove), y se cruza cada predicción con el
%   archivo de origen de la ventana (*_meta.csv).
%   -> resultados/analisis_confundido/*.csv

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
cfg = config_niveles();
weka_iniciar(cfg);
dirO = fullfile(cfg.rutas.resultados, 'analisis_confundido');
if ~isfolder(dirO), mkdir(dirO); end

% Escenarios: (ARFF, conjunto, modelo)
H = readtable(fullfile(cfg.rutas.resultados, 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), 'Delimiter', ',', 'TextType', 'char');
H.Properties.VariableNames = {'modelo','especificacion'};
H0 = readtable(fullfile(cfg.rutas.resultados, 'ajuste', 'hiperparametros_elegidos.csv'), 'Delimiter', ',', 'TextType', 'char');
H0.Properties.VariableNames = {'modelo','especificacion'};
mejor = strtrim(fileread(fullfile(cfg.rutas.resultados, 'calibracion', 'mejor_conjunto.txt')));
esc = {
  'Calibrado / MLP ajustado',   fullfile(cfg.rutas.resultados, 'arff', 'calibracion'), mejor, ...
        modelo_desde_id(H.especificacion{strcmp(H.modelo, 'MLP (ajustado)')});
  'Calibrado / SVM artículo',   fullfile(cfg.rutas.resultados, 'arff', 'calibracion'), mejor, ...
        especificar_modelo('SVM', struct());
  'Artículo HRS64+HO+HV / MLP ajustado', fullfile(cfg.rutas.resultados, 'arff', 'd512'), 'HRS64+HO+HV', ...
        modelo_desde_id(H0.especificacion{strcmp(H0.modelo, 'MLP (ajustado)')});
  'Artículo HRS64+HO+HV / SVM artículo', fullfile(cfg.rutas.resultados, 'arff', 'd512'), 'HRS64+HO+HV', ...
        especificar_modelo('SVM', struct())};
semillas = cfg.ajuste.semillas_final;
nw = usar_paralelo(cfg);

todas = table();
for e = 1:size(esc, 1)
    I = readtable(fullfile(esc{e,2}, 'indice.csv'), 'TextType', 'char');
    I = I(strcmp(I.condicion, 'KMI') & strcmp(I.esquema, '3c') & strcmp(I.alcance, 'intra'), :);
    partes = cell(height(I), 1);
    ruta_base = esc{e,2}; conj = esc{e,3}; m = esc{e,4};
    parfor (s = 1:height(I), nw)
        partes{s} = predicciones_sesion(fullfile(ruta_base, I.archivo{s}), conj, m, cfg, semillas);
    end
    P = vertcat(partes{:});
    P.escenario = repmat(esc(e,1), height(P), 1);
    todas = [todas; P]; %#ok<AGROW>
    fprintf('%s: %d predicciones\n', esc{e,1}, height(P));
end
writetable(todas, fullfile(dirO, 'predicciones.csv'));

%% Tablas: predicción del REPOSO según el archivo de origen
res = {};
for e = 1:size(esc, 1)
    P = todas(strcmp(todas.escenario, esc{e,1}), :);
    fprintf('\n=== %s ===\n', esc{e,1});
    for nivel = [10 40]
        R = P(P.clase == 0 & P.nivel_archivo == nivel, :);
        pct = 100 * [mean(strcmp(R.prediccion, 'N00')), mean(strcmp(R.prediccion, 'N10')), mean(strcmp(R.prediccion, 'N40'))];
        fprintf('  Reposo del archivo %2d %%  (n=%5d): -> N00 %5.1f %% | N10 %5.1f %% | N40 %5.1f %%\n', nivel, height(R), pct);
        res(end+1, :) = {esc{e,1}, nivel, height(R), pct(1), pct(2), pct(3)}; %#ok<AGROW>
    end
    % Entre los reposos clasificados como tarea: ¿coincide con su archivo?
    R = P(P.clase == 0 & ~strcmp(P.prediccion, 'N00'), :);
    coincide = (strcmp(R.prediccion, 'N10') & R.nivel_archivo == 10) | (strcmp(R.prediccion, 'N40') & R.nivel_archivo == 40);
    k = nnz(coincide); n = height(R);
    % Prueba exacta (binomial bilateral aproximada por normal; n grande) vs 50 %
    z = (k - n/2) / sqrt(n/4);
    p = erfc(abs(z)/sqrt(2));
    fprintf('  Reposos clasificados como tarea: %.1f %% coinciden con el nivel de su archivo (n=%d, p=%.2g; sin confundido ~50 %%)\n', 100*k/n, n, p);
    res(end+1, :) = {esc{e,1}, NaN, n, NaN, 100*k/n, p}; %#ok<AGROW>
end
T = cell2table(res, 'VariableNames', {'escenario','nivel_archivo_reposo','n','pct_N00','pct_N10_o_coincide','pct_N40_o_p'});
writetable(T, fullfile(dirO, 'resumen_reposo_por_archivo.csv'));
fprintf('\nResultados en %s\n', dirO);

%% ======================================================================
function P = predicciones_sesion(ruta, conjunto, m, cfg, semillas)
% CV 5-fold manual (idéntica a crossValidateModel) guardando cada predicción.
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
