function T = listar_registros(cfg)
%LISTAR_REGISTROS  Inventario de todos los .mat bajo cfg.rutas.datos.
%   T = listar_registros(cfg) devuelve una tabla con una fila por archivo:
%   ruta, archivo, paciente, mes, condicion (KMI/FIS), nivel (10/40),
%   formato ('A' o 'legado'), n_muestras, n_stims, estado y motivo.
%
%   Reglas (portadas de Niveles3/src/io/loader.py):
%   - paciente: carpeta 'Px.006' -> 'Px006'
%   - mes: primer entero del nombre de la carpeta ('1°MES' -> 1). NUNCA del
%     nombre de archivo (p. ej. px006_10_2.mat está en 3°MES).
%   - condición: 'img' en el nombre -> KMI (imaginería), si no FIS (físico)
%   - nivel: patrón _(10|40)(img)?_ del nombre; si difiere de max(stims) se
%     registra un aviso (stims no siempre se actualizó en la adquisición).
%   - formato A: EEG_all + 80 stims con patrón [0 0 0 0 L L L L]x10.
%     Formato legado (solo B): su temporización no es alineable -> excluido.

d = dir(fullfile(cfg.rutas.datos, '**', '*.mat'));
n = numel(d);
c = cell(n,1);
T = table(c, c, c, zeros(n,1), c, zeros(n,1), c, zeros(n,1), zeros(n,1), c, c, ...
    'VariableNames', {'ruta','archivo','paciente','mes','condicion','nivel', ...
    'formato','n_muestras','n_stims','estado','motivo'});

for k = 1:n
    ruta = fullfile(d(k).folder, d(k).name);
    % OJO: no usar fileparts con 'Px.006' (tomaría '.006' como extensión)
    partes = strsplit(d(k).folder, filesep);
    nombre_mes = partes{end};
    nombre_pac = partes{end-1};
    [~, base] = fileparts(ruta);

    T.ruta{k}     = ruta;
    T.archivo{k}  = base;
    T.paciente{k} = strrep(regexprep(nombre_pac, '(?i)^px', 'Px'), '.', '');
    tok = regexp(nombre_mes, '\d+', 'match', 'once');
    T.mes(k) = str2double(tok);
    if contains(lower(base), 'img'), T.condicion{k} = 'KMI'; else, T.condicion{k} = 'FIS'; end
    tok = regexp(lower(base), '_(10|40)(img)?_\d+$', 'tokens', 'once');
    if isempty(tok), T.nivel(k) = NaN; else, T.nivel(k) = str2double(tok{1}); end

    % Inspección del contenido (sin cargar todo cuando no hace falta)
    vars = whos('-file', ruta);
    nombres = {vars.name};
    S = load(ruta, 'stims');
    st = double(S.stims(:)');
    T.n_stims(k) = numel(st);
    if ismember('EEG_all', nombres)
        tam = vars(strcmp(nombres,'EEG_all')).size;
        T.formato{k} = 'A';
    else
        tam = vars(strcmp(nombres,'B')).size;
        T.formato{k} = 'legado';
    end
    T.n_muestras(k) = tam(2);

    [ok, motivo] = validar(T(k,:), st, tam, cfg);
    if ok
        T.estado{k} = 'ok'; T.motivo{k} = motivo;      % motivo = aviso (si lo hay)
    else
        T.estado{k} = 'excluido'; T.motivo{k} = motivo;
    end
end
T = sortrows(T, {'paciente','mes','condicion','nivel'});
end

function [ok, motivo] = validar(fila, st, tam, cfg)
ok = false; motivo = '';
if strcmp(fila.formato{1}, 'legado')
    if cfg.inventario.excluir_formato_legado
        motivo = 'formato legado (solo B): stims no alineable con la señal'; return
    end
    % Incluido solo para los análisis por etiqueta de archivo (10 vs 40):
    % sus stims no tienen temporización alineable (ver cargar_registro).
    if tam(1) ~= numel(cfg.canales.nombres) || isnan(fila.nivel)
        motivo = 'formato legado con canales o nivel no válidos'; return
    end
    ok = true; motivo = 'aviso: formato legado (solo B, ~89 s): sin bloques de stims, solo 10 vs 40';
    return
end
if tam(1) ~= numel(cfg.canales.nombres)
    motivo = sprintf('%d canales (se esperaban %d)', tam(1), numel(cfg.canales.nombres)); return
end
if numel(st) ~= cfg.inventario.n_stims_esperado
    motivo = sprintf('%d stims (se esperaban %d)', numel(st), cfg.inventario.n_stims_esperado); return
end
patron = repmat(cfg.inventario.patron_ciclo, 1, numel(st)/numel(cfg.inventario.patron_ciclo));
if ~isequal(st > 0, logical(patron))
    motivo = 'patrón de stims inesperado'; return
end
if isnan(fila.nivel)
    motivo = 'no se pudo leer el nivel (10/40) del nombre de archivo'; return
end
r = tam(2) / numel(st);
if r < cfg.inventario.muestras_por_stim(1) || r > cfg.inventario.muestras_por_stim(2)
    motivo = sprintf('N/K = %.1f fuera de rango', r); return
end
ok = true;
% El nivel se toma del NOMBRE (como Niveles3); stims solo da la temporización
% reposo/tarea. En algunas sesiones el programa de adquisición no actualizó
% el valor de stims (Px006-M3: todo 40; Px009-M3: todo 10) o lo intercambió
% (Px008-M3-KMI). Se deja constancia como aviso.
if max(st) ~= fila.nivel
    motivo = sprintf('aviso: nivel del nombre (%g) ~= max(stims) (%g)', fila.nivel, max(st));
end
end
