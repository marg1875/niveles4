%% CORRIDA_CONTROLES  Mismo pipeline de los pacientes en los 5 controles sanos (experimento).
%  Controles C01-C05 (Datos_controles/): mismo protocolo, mismo experimentador
%  y una sesión por persona. Se repite lo mínimo para las tablas de exactitud
%  con el SVM y sus matrices de confusión, con las mismas funciones y
%  parámetros que en los pacientes:
%  1) Etapa 1: inventario y preprocesamiento por sesión (filtro, ASR k = 30,
%     CAR, ICA), como en corrida_asr.
%  2) Tres clases dentro de la sesión (reposo, 10 y 40 %): HRS64 + HO + HV,
%     ventanas de 512 muestras, SVM del artículo, CV de 5 particiones con 10
%     semillas (etapas 2-3 y el escenario de fig_confusion).
%  3) Ventanas contiguas de 4 s (como corrida_legados) y las configuraciones
%     de la tabla final de los pacientes (validacion_final.csv) con el SVM
%     ajustado ig8: dentro de la sesión (etapa 17), LOSO con permutación
%     exacta (con una sesión por persona, la LOSO es también la LOPO), voto
%     por registro y decisión por sesión (etapa 16).
%  4) Matrices de confusión de pacientes y controles con las mismas funciones.
%  -> resultados_controles/

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_controles');
cfg = config_niveles();
assert(cfg.asr.activar && cfg.asr.k == 30 && ~cfg.artefactos.activar, 'config_niveles: se espera ASR k = 30 sin rechazo');
cfg.rutas.datos = cfg.rutas.datos_controles;
cfg.segmentacion.duraciones = 512;
cfg.experimentos.fase = 'A';
cfgP = cfg; cfgP.rutas.resultados = fullfile(raiz, 'resultados_asr');   % pacientes, para comparar
if ~isfolder(cfg.rutas.resultados), mkdir(cfg.rutas.resultados); end
save(fullfile(cfg.rutas.resultados, 'cfg_usada.mat'), 'cfg');
registrar(cfg, '==== Controles sanos: mismo pipeline que los pacientes ====');
nw = usar_paralelo(cfg);

%% 1) Preprocesamiento (una W por sesión)
etapa1_preprocesar(cfg);
ses = readtable(fullfile(cfg.rutas.resultados, 'sesiones.csv'), 'TextType', 'char', 'Delimiter', ',');

%% 2) Tres clases dentro de la sesión (HRS64 + HO + HV, SVM del artículo)
etapa2_caracteristicas(cfg);
etapa3_exportar_arff(cfg);
dirO = fullfile(cfg.rutas.resultados, 'analisis_confundido');
if ~isfolder(dirO), mkdir(dirO); end
esc = 'Artículo HRS64+HO+HV / SVM artículo';
archP = fullfile(dirO, 'predicciones.csv');
if ~isfile(archP)
    dirA = fullfile(cfg.rutas.resultados, 'arff', 'd512');
    I = readtable(fullfile(dirA, 'indice.csv'), 'TextType', 'char');
    I = I(strcmp(I.condicion, 'KMI') & strcmp(I.esquema, '3c') & strcmp(I.alcance, 'intra'), :);
    m3 = especificar_modelo('SVM', struct());
    partes = cell(height(I), 1); arch = I.archivo; sem = cfg.ajuste.semillas_final;
    parfor (s = 1:height(I), nw)
        partes{s} = predicciones_sesion(fullfile(dirA, arch{s}), 'HRS64+HO+HV', m3, cfg, sem);
    end
    P = vertcat(partes{:});
    P.escenario = repmat({esc}, height(P), 1);
    writetable(P, archP);
end
fuentes = {'pacientes', fullfile(cfgP.rutas.resultados, 'analisis_confundido', 'predicciones.csv'); 'controles', archP};
[A3, M3] = deal(table());
for g = 1:2
    P = readtable(fuentes{g, 2}, 'TextType', 'char');
    [a, s3, mu, de] = confusion3(P(strcmp(P.escenario, esc), :));
    A3 = [A3; table(repmat(fuentes(g, 1), numel(a), 1), s3, a, 'VariableNames', {'grupo','sesion','exactitud'})]; %#ok<AGROW>
    cl = {'reposo'; '10'; '40'};
    M3 = [M3; table(repmat(fuentes(g, 1), 3, 1), cl, mu(:, 1), de(:, 1), mu(:, 2), de(:, 2), mu(:, 3), de(:, 3), ...
        'VariableNames', {'grupo','real','pred_reposo','pred_reposo_de','pred_10','pred_10_de','pred_40','pred_40_de'})]; %#ok<AGROW>
    registrar(cfg, '  3c dentro de la sesión, %s: %.1f ± %.1f %% (%d sesiones) | diagonal reposo %.1f, 10 %.1f, 40 %.1f %%', ...
        fuentes{g, 1}, mean(a), std(a), numel(a), mu(1, 1), mu(2, 2), mu(3, 3));
end
writetable(A3, fullfile(dirO, 'exactitud_3c_por_sesion.csv'));
writetable(M3, fullfile(dirO, 'confusion_3c.csv'));

%% 3) Ventanas contiguas de 4 s (curvas de la etapa 13, como corrida_legados)
b = cfg.banda;
dirW = fullfile(cfg.rutas.resultados, b.dir_w);
if ~isfolder(dirW), mkdir(dirW); end
archW = fullfile(dirW, 'w_W1000.mat');
if ~isfile(archW)
    t0 = tic; L = 1000;
    Cand = candidatos_rango(512, cfg);
    baseH = struct('nombre', {}, 'modo', {}, 'm_max', {}, 'm_min', {});
    for r = {'Re64a2', 'Re128a4'}
        c = Cand(strcmp({Cand.nombre}, r{1}));
        baseH(end+1) = struct('nombre', c.nombre, 'modo', c.modo, 'm_max', c.m_max, 'm_min', c.m_min); %#ok<SAGROW>
    end
    [~, lay] = caract_registro(cfg, L, baseH);
    extractor = @(V) single(caract_curvas(V, cfg, lay));
    rech = struct('id', 'R0', 'umbral_z', Inf, 'min_canales', 1);
    k = logical(ses.completa) & strcmp(ses.condicion, 'KMI');
    sid = ses.sesion(k); cnd = ses.condicion(k); dirS = fullfile(cfg.rutas.resultados, 'senales');
    partes = cell(numel(sid), 1);
    parfor (s = 1:numel(sid), nw)
        partes{s} = ventanas_tiempo_sesion(fullfile(dirS, [sid{s} '.mat']), sid{s}, cnd{s}, baseH, L, cfg, ...
            rech, extractor, max(L, 512));
    end
    W = vertcat(partes{:});
    W.feat = single(matriz_caracteristicas(W));
    save(archW, 'W', '-v7.3');
    registrar(cfg, '  W1000: %d sesiones KMI, %d ventanas (%.1f min)', numel(sid), height(W), toc(t0)/60);
end

%% 3-4) Configuraciones finales: dentro de la sesión, LOSO, voto y matrices 2 x 2
dirB = fullfile(cfg.rutas.resultados, b.carpeta);
if ~isfolder(dirB), mkdir(dirB); end
Hc = readtable(fullfile(raiz, 'resultados', 'ajuste_calibrado', 'hiperparametros_elegidos.csv'), ...
    'Delimiter', ',', 'TextType', 'char');
Hc.Properties.VariableNames = {'modelo','especificacion'};
m = modelo_desde_id(Hc.especificacion{strcmp(Hc.modelo, 'SVM (ajustado)')});
Vp = readtable(fullfile(cfgP.rutas.resultados, b.carpeta, 'validacion_final.csv'), 'TextType', 'char');
Ip = readtable(fullfile(cfgP.rutas.resultados, b.carpeta, 'intra_sesion.csv'), 'TextType', 'char');
grupos = {'pacientes', cfgP; 'controles', cfg};
[Tab, Conf, Ses] = deal(table());
for g = 1:2
    [T, X, nombres, fam] = familias_banda(grupos{g, 2}, b);
    registrar(cfg, '  %s: %d ventanas KMI de 4 s, %d sesiones', grupos{g, 1}, height(T), numel(unique(T.sesion)));
    for c = 1:height(Vp)
        t0 = tic;
        j = ismember(fam, strsplit(Vp.conjunto{c}, '+'));
        Tk = T; Tk.feat = X(:, j); nom = nombres(j);
        [aI, KI] = intra_por_sesion(Tk, nom, m, cfg, 1:10, 5);
        pS = predicciones_grupo(Tk, nom, m, cfg, Tk.sesion);
        [aS, sesA] = acc_sesion(Tk, pS);
        [vA, vS] = votar(Tk, pS);
        if g == 1
            % pacientes: las cifras publicadas; aquí solo se comprueba que se reproducen
            registrar(cfg, '    %-36s intra %.1f (pub. %.1f) | LOSO %.1f (pub. %.1f)', Vp.conjunto{c}, ...
                mean(aI, 'omitnan'), Ip.intra(c), mean(aS), Vp.loso(c));
            ci = [Vp.loso_ic_inf(c) Vp.loso_ic_sup(c)]; p = Vp.p_perm_sesion_1000(c);
        else
            % controles: permutación exacta de las 2^5 combinaciones de sesiones volteadas
            sesC = unique(Tk.sesion, 'stable'); nS = numel(sesC); nula = zeros(2^nS, 1);
            parfor (q = 1:2^nS, nw)
                vol = sesC(bitget(q - 1, 1:nS) == 1);
                nula(q) = mean(acc_sesion(voltear(Tk, vol), predicciones_grupo(Tk, nom, m, cfg, Tk.sesion, vol)));
            end
            p = mean(nula >= mean(aS));
            ci = boot(aS, regexprep(sesA, '_.*$', ''));
            registrar(cfg, '    %-36s intra %.1f %% | LOSO %.1f %% [%.0f-%.0f] p = %.3f | voto %.0f %% | sesión %.0f %% (%.1f min)', ...
                Vp.conjunto{c}, mean(aI, 'omitnan'), mean(aS), ci, p, vA(1), vS(1), toc(t0)/60);
        end
        Tab = [Tab; table(grupos(g, 1), Vp.conjunto(c), nnz(j), mean(aI, 'omitnan'), std(aI, 'omitnan'), ...
            mean(aS), ci(1), ci(2), nnz(aS > 50), numel(aS), p, vA(1), vS(1), vS(2), ...
            'VariableNames', {'grupo','conjunto','n_atributos','intra','intra_de','loso','loso_ic_inf','loso_ic_sup', ...
            'sobre50','n_ses','p','voto_archivo','decision_sesion','decision_sesion_p'})]; %#ok<AGROW>
        % matrices 2 x 2 (% de cada clase real, media y DE entre sesiones)
        KL = zeros(2, 2, numel(sesA));
        for s = 1:numel(sesA)
            k = strcmp(Tk.sesion, sesA{s}) & isfinite(pS);
            KL(:, :, s) = accumarray([(Tk.nivel(k) == 40) + 1, (pS(k) == 40) + 1], 1, [2 2]);
        end
        for v = {'intra', KI(:, :, isfinite(aI)); 'LOSO', KL}'
            Pc = 100 * v{2} ./ sum(v{2}, 2);
            mu = mean(Pc, 3); de = std(Pc, 0, 3);
            Conf = [Conf; table(repmat(grupos(g, 1), 2, 1), repmat(Vp.conjunto(c), 2, 1), repmat(v(1), 2, 1), [10; 40], ...
                mu(:, 1), de(:, 1), mu(:, 2), de(:, 2), 'VariableNames', ...
                {'grupo','conjunto','validacion','real','pred_10','pred_10_de','pred_40','pred_40_de'})]; %#ok<AGROW>
        end
        Ses = [Ses; table(repmat(grupos(g, 1), numel(sesA), 1), repmat(Vp.conjunto(c), numel(sesA), 1), sesA, aS, ...
            aI(ismember(unique(Tk.sesion, 'stable'), sesA)), 'VariableNames', {'grupo','conjunto','sesion','loso','intra'})]; %#ok<AGROW>
    end
end
writetable(Tab, fullfile(dirB, 'tabla_pacientes_controles.csv'));
writetable(Conf, fullfile(dirB, 'confusion_2c.csv'));
writetable(Ses, fullfile(dirB, 'por_sesion.csv'));
registrar(cfg, '==== Controles listos: %s ====', cfg.rutas.resultados);

%% ======================================================================
function [acc, ses, mu, de] = confusion3(P)
% Exactitud por sesión y matriz 3 x 3 (% de cada clase real, media y DE entre sesiones).
cl = {'N00', 'N10', 'N40'};
ses = unique(P.sesion);
C = nan(numel(ses), 3, 3); acc = nan(numel(ses), 1);
for s = 1:numel(ses)
    i = strcmp(P.sesion, ses{s});
    rc = (P.clase(i) == 10) + 2 * (P.clase(i) == 40) + 1;
    pc = cellfun(@(x) find(strcmp(cl, x)), P.prediccion(i));
    K = accumarray([rc pc], 1, [3 3]);
    C(s, :, :) = 100 * K ./ sum(K, 2);
    acc(s) = 100 * mean(rc == pc);
end
mu = squeeze(mean(C, 1)); de = squeeze(std(C, 0, 1));
end

function T = voltear(T, vol)
v = ismember(T.sesion, vol); T.nivel(v) = 50 - T.nivel(v);
end

function [a, ses] = acc_sesion(T, pred)
ses = unique(T.sesion, 'stable'); a = nan(numel(ses), 1);
for s = 1:numel(ses)
    k = strcmp(T.sesion, ses{s}) & isfinite(pred);
    if numel(unique(T.nivel(k))) == 2, a(s) = 100 * mean(pred(k) == T.nivel(k)); end
end
ok = isfinite(a); a = a(ok); ses = ses(ok);
end

function ci = boot(a, pac)
rs = RandStream('mt19937ar', 'Seed', 7); up = unique(pac); mm = zeros(10000, 1);
for i = 1:10000
    sel = up(randi(rs, numel(up), numel(up), 1));
    mm(i) = mean(cell2mat(cellfun(@(p) a(strcmp(pac, p)), sel, 'UniformOutput', false)));
end
ci = quantile(mm, [0.025 0.975]);
end

function [A, S] = votar(T, pred)
[g, ses, ~] = findgroups(T.sesion, T.archivo);
f40 = splitapply(@(p) mean(p(isfinite(p)) == 40), pred, g);
niv = splitapply(@(x) x(1), T.nivel, g);
ac = (f40 > 0.5) == (niv == 40);
A = [100 * mean(ac), 1 - binocdf(nnz(ac) - 1, numel(ac), 0.5)];
us = unique(ses, 'stable'); d = nan(numel(us), 1);
for s = 1:numel(us)
    i10 = strcmp(ses, us{s}) & niv == 10; i40 = strcmp(ses, us{s}) & niv == 40;
    if any(i10) && any(i40), d(s) = f40(i40) > f40(i10); end
end
d = d(isfinite(d));
S = [100 * mean(d), 1 - binocdf(nnz(d) - 1, numel(d), 0.5)];
end
