function etapa5_reportes(cfg)
%ETAPA5_REPORTES  Tablas y figuras a partir de resultados/weka/metricas.csv.
%   Fase A: precisión vs duración (configuración de referencia).
%   Fase B: precisión por configuración de canales (cajas por sesión) +
%           Friedman (sesiones = bloques) + Wilcoxon por pares con Holm.
%   Análogos del artículo en la configuración y duración de referencia:
%     Tabla 1 (por partición p), Tabla 2 (HRS/HO/HV/combinado), Tabla 3
%     (TPR, FPR, precisión, F, ROC), Fig. 7 (HV por nivel), Fig. 8
%     (precisión por sesión y matriz de confusión).
%   Extra: precisión por mes (longitudinal) y chequeo de fuga.
%   Las medias "intra" promedian las sesiones (equivalente al promedio por
%   sujeto del artículo). El "mejor de la rejilla" es optimista: se elige
%   sobre los mismos folds con que se evalúa.
%   -> resultados/reportes/*.csv y resultados/reportes/figuras/*.png

dirW = fullfile(cfg.rutas.resultados, 'weka');
dirR = fullfile(cfg.rutas.resultados, 'reportes');
dirF = fullfile(dirR, 'figuras');
if ~isfolder(dirF), mkdir(dirF); end
R = readtable(fullfile(dirW, 'metricas.csv'), 'TextType', 'char');
R = R(strcmp(R.estado, 'ok'), :);
E = R(strcmp(R.cv_modo, 'estratificada'), :);
ref = cfg.canales.config_referencia;
dref = cfg.segmentacion.duracion_referencia;
crit = cfg.experimentos.criterio;
mejor = dref;
if isfile(fullfile(dirW, 'mejor_duracion.txt')), mejor = load(fullfile(dirW, 'mejor_duracion.txt')); end
clfs = unique(E.clasificador, 'stable');

%% ---------------- Fase A: precisión vs duración
A = E(strcmp(E.config, ref), :);
TA = resumir(A, {'condicion','esquema','alcance','conjunto','clasificador','duracion'});
TA.segundos = TA.duracion / cfg.fs;
writetable(TA, fullfile(dirR, 'faseA_precision_vs_duracion.csv'));
for cond = cfg.experimentos.condiciones
    for esq = cfg.experimentos.esquemas
        f = figure('Visible', 'off', 'Position', [100 100 1100 400]);
        conjs = {'HRS64+HO+HV', 'HRS64', 'HO', 'HV'};
        for ic = 1:numel(conjs)
            subplot(1, numel(conjs), ic); hold on
            for c = 1:numel(clfs)
                s = TA(strcmp(TA.condicion,cond{1}) & strcmp(TA.esquema,esq{1}) & ...
                    strcmp(TA.alcance,'intra') & strcmp(TA.conjunto,conjs{ic}) & ...
                    strcmp(TA.clasificador,clfs{c}), :);
                s = sortrows(s, 'segundos');
                errorbar(s.segundos, s.acc_media, s.acc_de, '-o', 'DisplayName', clfs{c});
            end
            set(gca, 'XScale', 'log'); grid on
            xlabel('Duración de ventana (s)'); ylabel('Precisión intra-sujeto (%)');
            title(conjs{ic});
            if ic == 1, legend('Location', 'best'); end
        end
        sgtitle(sprintf('Fase A - %s %s (config %s): media \\pm DE entre sesiones', cond{1}, esq{1}, ref));
        exportgraphics(f, fullfile(dirF, sprintf('faseA_%s_%s.png', cond{1}, esq{1})), 'Resolution', 150);
        close(f);
    end
end

%% ---------------- Fase B: configuraciones de canales
estad = {};
for d = unique([dref, mejor])
    B = E(E.duracion == d & strcmp(E.condicion, crit.condicion) & strcmp(E.esquema, crit.esquema) & ...
        strcmp(E.alcance, 'intra') & strcmp(E.conjunto, crit.conjunto), :);
    cfgs = intersect({cfg.canales.configs.id}, unique(B.config), 'stable');
    if numel(cfgs) < 2, continue, end
    TB = resumir(E(E.duracion == d, :), {'condicion','esquema','alcance','conjunto','clasificador','config'});
    writetable(TB, fullfile(dirR, sprintf('faseB_configs_d%d.csv', d)));
    f = figure('Visible', 'off', 'Position', [100 100 1300 420]);
    for c = 1:numel(clfs)
        Bc = B(strcmp(B.clasificador, clfs{c}), :);
        ses = unique(Bc.id);
        M = nan(numel(ses), numel(cfgs));
        for j = 1:numel(cfgs)
            for s = 1:numel(ses)
                v = Bc.acc(strcmp(Bc.config, cfgs{j}) & strcmp(Bc.id, ses{s}));
                if ~isempty(v), M(s,j) = v(1); end
            end
        end
        [p, chi2] = prueba_friedman(M);
        estad(end+1,:) = {d, clfs{c}, size(M,1), chi2, p}; %#ok<AGROW>
        W = wilcoxon_holm(M, cfgs);
        writetable(W, fullfile(dirR, sprintf('faseB_wilcoxon_d%d_%s.csv', d, clfs{c})));
        subplot(1, numel(clfs), c);
        grupo = categorical(repmat(cfgs, size(M,1), 1), cfgs);
        boxchart(grupo(:), M(:)); grid on
        ylabel('Precisión intra-sujeto (%)');
        title(sprintf('%s (Friedman p = %.3g)', clfs{c}, p));
    end
    sgtitle(sprintf('Fase B - %s %s, %s, %d muestras (%.2f s)', crit.condicion, crit.esquema, ...
        crit.conjunto, d, d/cfg.fs));
    exportgraphics(f, fullfile(dirF, sprintf('faseB_configs_d%d.png', d)), 'Resolution', 150);
    close(f);
    % Curva precisión vs número de canales (configuraciones escalonadas E1..E9)
    esc = cfgs(startsWith(cfgs, 'E'));
    nch = arrayfun(@(i) numel(cfg.canales.configs(strcmp({cfg.canales.configs.id}, esc{i})).canales), 1:numel(esc));
    f = figure('Visible', 'off'); hold on
    for c = 1:numel(clfs)
        m = arrayfun(@(i) mean(B.acc(strcmp(B.config, esc{i}) & strcmp(B.clasificador, clfs{c}))), 1:numel(esc));
        plot(nch, m, '-o', 'DisplayName', clfs{c});
    end
    grid on; xlabel('Número de canales (E1 \rightarrow E9)'); ylabel('Precisión intra-sujeto media (%)');
    legend('Location', 'best'); title(sprintf('Configuraciones escalonadas, %d muestras', d));
    exportgraphics(f, fullfile(dirF, sprintf('faseB_escalonadas_d%d.png', d)), 'Resolution', 150);
    close(f);
end
if ~isempty(estad)
    writetable(cell2table(estad, 'VariableNames', {'duracion','clasificador','n_sesiones','chi2','p_friedman'}), ...
        fullfile(dirR, 'faseB_friedman.csv'));
end

%% ---------------- Tablas 1-3 análogas (config y duración de referencia)
Rf = E(strcmp(E.config, ref) & E.duracion == dref, :);
for cond = cfg.experimentos.condiciones
    for esq = cfg.experimentos.esquemas
        S = Rf(strcmp(Rf.condicion, cond{1}) & strcmp(Rf.esquema, esq{1}), :);
        if isempty(S), continue, end
        sufijo = sprintf('%s_%s', cond{1}, esq{1});
        % Tabla 1: HRSp+HO+HV por p
        t1 = resumir(S(startsWith(S.conjunto, 'HRS') & endsWith(S.conjunto, '+HO+HV'), :), ...
            {'alcance','conjunto','clasificador'});
        writetable(t1, fullfile(dirR, sprintf('tabla1_particiones_%s.csv', sufijo)));
        % Tabla 2: métodos
        t2 = resumir(S(ismember(S.conjunto, {'HRS64','HO','HV','HRS64+HO+HV'}), :), ...
            {'alcance','conjunto','clasificador'});
        writetable(t2, fullfile(dirR, sprintf('tabla2_metodos_%s.csv', sufijo)));
        % Tabla 3: métricas ponderadas (HRS64+HO+HV, intra: media +- DE)
        s3 = S(strcmp(S.conjunto, 'HRS64+HO+HV') & strcmp(S.alcance, 'intra'), :);
        [g, cl] = findgroups(s3.clasificador);
        f = @(x) mean(x, 'omitnan'); sd = @(x) std(x, 'omitnan');
        t3 = table(cl, 100*splitapply(f, s3.wtpr, g), 100*splitapply(sd, s3.wtpr, g), ...
            100*splitapply(f, s3.wfpr, g), 100*splitapply(sd, s3.wfpr, g), ...
            100*splitapply(f, s3.wprec, g), 100*splitapply(sd, s3.wprec, g), ...
            100*splitapply(f, s3.wf, g), 100*splitapply(sd, s3.wf, g), ...
            100*splitapply(f, s3.wroc, g), 100*splitapply(sd, s3.wroc, g), ...
            'VariableNames', {'clasificador','TPR','TPR_de','FPR','FPR_de','precision','precision_de', ...
            'F','F_de','ROC','ROC_de'});
        writetable(t3, fullfile(dirR, sprintf('tabla3_metricas_%s.csv', sufijo)));
        % Fig. 8 análoga: precisión por sesión (kNN) + matriz de confusión media
        k8 = s3(strcmp(s3.clasificador, 'kNN'), :);
        if ~isempty(k8)
            fg = figure('Visible', 'off', 'Position', [100 100 1000 380]);
            subplot(1,2,1); bar(categorical(k8.id), k8.acc); ylabel('Precisión (%)'); grid on
            title(sprintf('kNN por sesión - %s', sufijo), 'Interpreter', 'none');
            C = zeros(size(eval(k8.confusion{1})));
            for i = 1:height(k8)
                Ci = eval(k8.confusion{i});
                C = C + Ci ./ sum(Ci, 2);
            end
            C = 100 * C / height(k8);
            subplot(1,2,2);
            et = {'N00','N10','N40'}; if size(C,1) == 2, et = {'N10','N40'}; end
            heatmap(et, et, round(C, 1), 'XLabel', 'Predicha', 'YLabel', 'Real', ...
                'Title', 'Confusión media (% por fila)');
            exportgraphics(fg, fullfile(dirF, sprintf('fig8_%s.png', sufijo)), 'Resolution', 150);
            close(fg);
        end
    end
end

%% ---------------- Fig. 7 análoga: HV por nivel
archC = fullfile(cfg.rutas.resultados, 'caracteristicas', sprintf('caracteristicas_d%d.csv', dref));
if isfile(archC)
    C = readtable(archC, 'TextType', 'char');
    C = C(logical(C.sesion_completa), :);
    hv = startsWith(C.Properties.VariableNames, 'HV_');
    for cond = cfg.experimentos.condiciones
        S = C(strcmp(C.condicion, cond{1}), :);
        if isempty(S), continue, end
        fg = figure('Visible', 'off', 'Position', [100 100 1100 400]);
        subplot(1,2,1); hold on
        niveles = [0 10 40];
        for n = 1:3
            m = mean(S{S.clase == niveles(n), hv}, 1, 'omitnan');
            plot(1:nnz(hv), m, 'o', 'DisplayName', sprintf('%d %%', niveles(n)));
        end
        set(gca, 'XTick', 1:nnz(hv), 'XTickLabel', cfg.canales.nombres); grid on
        ylabel('HV medio'); legend; title(sprintf('HV por electrodo y nivel - %s', cond{1}));
        subplot(1,2,2);
        boxchart(categorical(S.clase), mean(S{:, hv}, 2, 'omitnan')); grid on
        xlabel('Nivel (% FVM)'); ylabel('HV (media de canales por ventana)');
        title('HV por nivel');
        exportgraphics(fg, fullfile(dirF, sprintf('fig7_hv_%s.png', cond{1})), 'Resolution', 150);
        close(fg);
        resumen = groupsummary(S(:, [{'clase'}, C.Properties.VariableNames(hv)]), 'clase', {'mean','std'});
        writetable(resumen, fullfile(dirR, sprintf('hv_por_nivel_%s.csv', cond{1})));
    end
end

%% ---------------- Longitudinal por mes y chequeo de fuga
L = E(strcmp(E.alcance, 'intra') & strcmp(E.config, ref) & E.duracion == dref, :);
if ~isempty(L)
    L.mes = cellfun(@(s) sscanf(regexp(s, 'M\d+', 'match', 'once'), 'M%d'), L.id);
    TL = resumir(L, {'condicion','esquema','conjunto','clasificador','mes'});
    writetable(TL, fullfile(dirR, 'longitudinal_por_mes.csv'));
end
G = R(strcmp(R.cv_modo, 'agrupada'), :);
if ~isempty(G)
    claves = {'duracion','dataset','config','conjunto','clasificador'};
    Eg = innerjoin(G(:, [claves, {'acc'}]), E(:, [claves, {'acc'}]), 'Keys', claves);
    Eg.Properties.VariableNames(end-1:end) = {'acc_agrupada','acc_estratificada'};
    Eg.diferencia = Eg.acc_estratificada - Eg.acc_agrupada;
    writetable(Eg, fullfile(dirR, 'chequeo_fuga.csv'));
    registrar(cfg, 'Fuga: CV estratificada supera a la agrupada en %.1f puntos en promedio.', mean(Eg.diferencia));
end
registrar(cfg, 'Etapa 5: reportes en %s', dirR);
end

function T = resumir(R, claves)
% Media, DE y n de las métricas principales agrupando por las claves.
[g, T] = findgroups(R(:, claves));
T.n = splitapply(@numel, R.acc, g);
T.acc_media = splitapply(@(x) mean(x, 'omitnan'), R.acc, g);
T.acc_de = splitapply(@(x) std(x, 'omitnan'), R.acc, g);
T.kappa_media = splitapply(@(x) mean(x, 'omitnan'), R.kappa, g);
T.mae_media = splitapply(@(x) mean(x, 'omitnan'), R.mae, g);
end
