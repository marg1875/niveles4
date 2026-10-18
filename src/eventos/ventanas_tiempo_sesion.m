function T = ventanas_tiempo_sesion(rutaS, sid, cond, fams, L, cfg, rechazos, extractor, paso)
%VENTANAS_TIEMPO_SESION  Ventanas contiguas de todo el registro con su tiempo y validez.
%   T = ventanas_tiempo_sesion(rutaS, sid, cond, fams, L, cfg, rechazos, extractor, paso)
%   Ventanas no traslapadas de L muestras (desde 0.5 s) de cada archivo de
%   la sesión, con: nivel (corregido: si los dos archivos tienen stims
%   distintos 10/40 se usa stims), minuto dentro del archivo y de la sesión
%   (marca "Created on" del .mat), características (caracterizadores_limites)
%   y una columna lógica 'valida' por variante de rechazo (rechazos(k):
%   struct con umbral_z y min_canales; umbral_z = Inf -> sin rechazo).
%   extractor (opcional): función V (C x L) -> fila de características; por
%   defecto caracterizadores_limites(V, fams, cfg).
%   paso (opcional): separación entre inicios de ventana (por defecto L,
%   ventanas contiguas). Con paso > L quedan huecos; con paso < L, traslape.
if nargin < 8 || isempty(extractor), extractor = @(V) caracterizadores_limites(V, fams, cfg); end
if nargin < 9, paso = L; end
S = load(rutaS);
fs = cfg.fs;
nivel = S.niveles;
stimsMax = cellfun(@max, S.stims);
if numel(S.X) == 2 && all(sort(stimsMax) == [10 40]), nivel = stimsMax; end
hora = NaT(1, numel(S.X));
for j = 1:numel(S.X), hora(j) = hora_mat(S.archivos{j}); end
dur = cellfun(@(x) size(x, 2), S.X) / fs;
t0 = min(hora - seconds(dur));
filas = {};
for j = 1:numel(S.X)
    X = double(S.X{j});
    N = size(X, 2);
    off = minutes(hora(j) - seconds(N / fs) - t0);
    [~, arch] = fileparts(S.archivos{j});
    inicios = round(0.5 * fs) + 1 : paso : N - L + 1;
    val = true(numel(inicios), numel(rechazos));
    for r = 1:numel(rechazos)
        if isinf(rechazos(r).umbral_z), continue, end
        c = cfg; c.artefactos.activar = true;
        c.artefactos.umbral_z = rechazos(r).umbral_z;
        c.artefactos.min_canales = rechazos(r).min_canales;
        val(:, r) = ventana_valida(marcar_artefactos(X, c), inicios, L).';
    end
    for i = 1:numel(inicios)
        s = inicios(i);
        f = extractor(X(:, s:s+L-1));
        filas(end+1, :) = {sid, cond, arch, nivel(j), (s + L/2) / fs / 60, off + (s + L/2) / fs / 60, ...
            f, val(i, :)}; %#ok<AGROW>
    end
end
T = cell2table(filas, 'VariableNames', {'sesion','condicion','archivo','nivel','t_archivo_min','t_min','feat','valida'});
end

function h = hora_mat(ruta)
fid = fopen(ruta, 'r'); b = fread(fid, 116, '*char')'; fclose(fid);
tok = regexp(b, 'Created on: (\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})', 'tokens', 'once');
h = datetime(regexprep(tok{1}, ' +', ' '), 'InputFormat', 'eee MMM d HH:mm:ss yyyy', 'Locale', 'en_US');
end
