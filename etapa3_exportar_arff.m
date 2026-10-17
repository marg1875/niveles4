function etapa3_exportar_arff(cfg)
%ETAPA3_EXPORTAR_ARFF  Construye los datasets ARFF para Weka.
%   Para cada duración L se generan los datasets:
%     condición {KMI, FIS} x esquema {3c: N00/N10/N40, 2c: N10/N40}
%     x alcance {intra: una sesión (equivale al "por sujeto" del artículo),
%                mes: todas las sesiones de un mes, global: todas}
%   Solo entran sesiones COMPLETAS (archivo de 10 % y de 40 %).
%   La clase N00 sale de los bloques de reposo de ambos archivos de la
%   sesión, submuestreados de forma determinista (equiespaciados) para
%   igualar el número de ventanas de N10/N40 (cfg.segmentacion.balancear_reposo).
%   Cada ARFF lleva los 16 canales x 9 características; las configuraciones
%   de canales y los conjuntos de características se seleccionan en Weka.
%   -> resultados/arff/d<L>/{intra,mes,global}/<COND>_<esq>_<id>.arff
%      + <...>_meta.csv (metadatos alineados por fila) + d<L>/indice.csv

dirC = fullfile(cfg.rutas.resultados, 'caracteristicas');
nombres = nombres_caracteristicas(cfg);
for L = cfg.segmentacion.duraciones
    dirA = fullfile(cfg.rutas.resultados, 'arff', sprintf('d%d', L));
    if isfile(fullfile(dirA, 'indice.csv'))
        registrar(cfg, 'Etapa 3: ARFF d=%d ya existen.', L); continue
    end
    archivoC = fullfile(dirC, sprintf('caracteristicas_d%d.csv', L));
    if ~isfile(archivoC), continue, end
    T = readtable(archivoC, 'TextType', 'char');
    T = T(logical(T.sesion_completa), :);
    T = balancear_reposo(T, cfg);
    T.etiqueta = arrayfun(@(c) sprintf('N%02d', c), T.clase, 'UniformOutput', false);

    indice = {};
    for ic = 1:numel(cfg.experimentos.condiciones)
        cond = cfg.experimentos.condiciones{ic};
        Tc = T(strcmp(T.condicion, cond), :);
        for ie = 1:numel(cfg.experimentos.esquemas)
            esq = cfg.experimentos.esquemas{ie};
            if strcmp(esq, '3c')
                clases = {'N00','N10','N40'}; Te = Tc;
            else
                clases = {'N10','N40'}; Te = Tc(Tc.clase > 0, :);
            end
            for ia = 1:numel(cfg.experimentos.alcances)
                alc = cfg.experimentos.alcances{ia};
                switch alc
                    case 'intra',  [grp, ids] = findgroups(Te.sesion);
                    case 'mes',    [grp, ms] = findgroups(Te.mes); ids = arrayfun(@(m) sprintf('M%d', m), ms, 'UniformOutput', false);
                    case 'global', grp = ones(height(Te),1); ids = {'GLOBAL'};
                end
                for k = 1:numel(ids)
                    D = Te(grp == k, :);
                    if isempty(D), continue, end
                    nombre = sprintf('%s_%s_%s', cond, esq, regexprep(ids{k}, ['_' cond '$'], ''));
                    ruta = fullfile(dirA, alc, [nombre '.arff']);
                    com = {sprintf('niveles4 (%s) - %s', cfg.version, datestr(now, 'yyyy-mm-dd HH:MM')), ... %#ok<TNOW1,DATST>
                        sprintf('Dataset %s | duracion %d muestras (%.3f s) | alcance %s', nombre, L, L/cfg.fs, alc), ...
                        sprintf('Filtro %s %g-%g Hz + notch %g Hz | ICA por sesion (sensores limpios)', ...
                            cfg.filtro.tipo, cfg.filtro.banda, cfg.filtro.notch_hz), ...
                        sprintf('Hurst: HO (bloque>=%d), HRS p=%s, HV (R/S del semivariograma, p=%d), <R/S> %s', ...
                            cfg.hurst.ho_bloque_min, mat2str(cfg.hurst.particiones_hrs), cfg.hurst.hv_particiones, cfg.hurst.promedio_rs), ...
                        'Canales: Fp1 Fp3 F3 Fz F4 Cz Pz P5 P3 O1 Oz T8 P8 P6 P4 T7 (g.Nautilus, 250 Hz)'};
                    exportar_arff(ruta, nombre, nombres, D{:, nombres}, clases, D.etiqueta, com);
                    meta = D(:, {'sesion','paciente','mes','condicion','archivo','nivel_archivo','clase','bloque','ancla'});
                    meta.grupo = strcat(meta.archivo, '_b', string(meta.bloque));   % para CV agrupada
                    writetable(meta, strrep(ruta, '.arff', '_meta.csv'), 'Encoding', 'UTF-8');
                    cuenta = arrayfun(@(c) nnz(strcmp(D.etiqueta, c)), clases);
                    indice(end+1, :) = {nombre, cond, esq, alc, ids{k}, height(D), mat2str(cuenta), ...
                        fullfile(alc, [nombre '.arff'])}; %#ok<AGROW>
                end
            end
        end
    end
    I = cell2table(indice, 'VariableNames', {'dataset','condicion','esquema','alcance','id', ...
        'n_instancias','n_por_clase','archivo'});
    writetable(I, fullfile(dirA, 'indice.csv'), 'Encoding', 'UTF-8');
    registrar(cfg, 'Etapa 3: d=%d -> %d datasets ARFF.', L, height(I));
end
end

function T = balancear_reposo(T, cfg)
% Submuestreo determinista de las ventanas de reposo por sesión: se eligen
% índices equiespaciados (ordenadas por archivo, bloque y ancla) hasta igualar
% el promedio de ventanas de N10 y N40 de la sesión.
if ~cfg.segmentacion.balancear_reposo, return, end
mantener = true(height(T), 1);
[g, ~] = findgroups(T.sesion);
for s = 1:max(g)
    idx = find(g == s);
    reposo = idx(T.clase(idx) == 0);
    objetivo = round(mean([nnz(T.clase(idx) == 10), nnz(T.clase(idx) == 40)]));
    if numel(reposo) > objetivo
        [~, o] = sortrows(T(reposo, {'archivo','bloque','ancla'}));
        reposo = reposo(o);
        elegidos = reposo(unique(round(linspace(1, numel(reposo), objetivo))));
        mantener(setdiff(reposo, elegidos)) = false;
    end
end
T = T(mantener, :);
end
