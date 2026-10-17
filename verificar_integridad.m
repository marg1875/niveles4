%% VERIFICAR_INTEGRIDAD  Revisión de integridad de todos los .mat de Datos/.
%
%   Para cada archivo: variables y tamaños, duración, NaN/Inf, canales planos,
%   saturación (valores repetidos), amplitud extrema, canales duplicados
%   (|r| > 0.999), potencia de línea (60 Hz) relativa, relación de B con
%   EEG_all, y contenido de stims/resps/respsK. Los criterios son
%   descriptivos: se marcan banderas para revisar, no se excluye nada.
%   -> resultados/integridad/integridad_archivos.csv, integridad_canales.csv

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
cfg = config_niveles();
dirI = fullfile(cfg.rutas.resultados, 'integridad');
if ~isfolder(dirI), mkdir(dirI); end
d = dir(fullfile(cfg.rutas.datos, '**', '*.mat'));
fs = cfg.fs;
F = {}; Cn = {};
for k = 1:numel(d)
    ruta = fullfile(d(k).folder, d(k).name);
    partes = strsplit(d(k).folder, filesep);
    v = whos('-file', ruta);
    S = load(ruta);
    vars = strjoin(arrayfun(@(x) sprintf('%s[%s]%s', x.name, strjoin(string(x.size), 'x'), x.class), v, 'UniformOutput', false), ' ');
    if isfield(S, 'EEG_all'), X = double(S.EEG_all); fuente = 'EEG_all'; else, X = double(S.B); fuente = 'B'; end
    [nc, N] = size(X);
    nNaN = nnz(isnan(X)); nInf = nnz(isinf(X));
    X(~isfinite(X)) = 0;
    sd = std(X, 0, 2);
    plano = sd < 1e-6;
    % saturación: fracción de muestras iguales a la anterior (valores pegados)
    pegado = mean(diff(X, 1, 2) == 0, 2);
    rango = max(X, [], 2) - min(X, [], 2);
    R = corrcoef(X.'); R(1:nc+1:end) = 0;
    dup = any(abs(R) > 0.999, 2);
    % potencia de línea relativa (59-61 Hz vs 1-45 Hz)
    P = abs(fft(X - mean(X, 2), [], 2)).^2; fr = (0:N-1) * fs / N;
    linea = sum(P(:, fr >= 59 & fr <= 61), 2) ./ sum(P(:, fr >= 1 & fr <= 45), 2);
    % B como cola de EEG_all
    bcola = NaN;
    if isfield(S, 'EEG_all') && isfield(S, 'B') && size(S.B, 2) <= N
        bcola = max(abs(S.EEG_all(:, end-size(S.B,2)+1:end) - S.B), [], 'all');
    end
    st = double(S.stims(:)');
    ex = '';
    if isfield(S, 'resps'), ex = sprintf('resps: %d valores, %d distintos de 0', numel(S.resps), nnz(S.resps)); end
    if isfield(S, 'respsK'), ex = [ex sprintf('; respsK: %s', strjoin(unique(string(S.respsK)), ','))]; end
    banderas = {};
    if nNaN || nInf, banderas{end+1} = 'NaN/Inf'; end
    if any(plano), banderas{end+1} = sprintf('%d canal(es) plano(s)', nnz(plano)); end
    if any(pegado > 0.05), banderas{end+1} = sprintf('%d canal(es) con >5%% de muestras repetidas', nnz(pegado > 0.05)); end
    if any(dup), banderas{end+1} = sprintf('%d canal(es) duplicado(s)', nnz(dup)); end
    if any(linea > 0.5), banderas{end+1} = sprintf('%d canal(es) con 60 Hz dominante', nnz(linea > 0.5)); end
    if any(rango > 10 * median(rango)), banderas{end+1} = sprintf('%d canal(es) con amplitud >10x la mediana', nnz(rango > 10 * median(rango))); end
    F(end+1, :) = {partes{end-1}, partes{end}, d(k).name, vars, fuente, nc, N, N/fs, nNaN + nInf, ...
        numel(st), mat2str(unique(st)), bcola, ex, strjoin(banderas, '; ')}; %#ok<AGROW>
    for c = 1:nc
        Cn(end+1, :) = {d(k).name, cfg.canales.nombres{min(c, end)}, sd(c), rango(c), pegado(c), linea(c), dup(c)}; %#ok<AGROW>
    end
end
T = cell2table(F, 'VariableNames', {'paciente','mes','archivo','variables','fuente_eeg','canales', ...
    'muestras','duracion_s','nan_inf','n_stims','valores_stims','max_dif_B_vs_cola','resps','banderas'});
writetable(T, fullfile(dirI, 'integridad_archivos.csv'));
writetable(cell2table(Cn, 'VariableNames', {'archivo','canal','desv_std','rango','frac_repetidas', ...
    'potencia_60Hz_rel','duplicado'}), fullfile(dirI, 'integridad_canales.csv'));
disp(T(:, {'paciente','mes','archivo','fuente_eeg','duracion_s','n_stims','valores_stims','max_dif_B_vs_cola','banderas'}));
fprintf('\nEsquemas de variables distintos:\n');
esq = regexprep(T.variables, '\[[^\]]*\]', '');
[u, ~, g] = unique(esq);
for i = 1:numel(u), fprintf('  %3d archivos: %s\n', nnz(g == i), u{i}); end
fprintf('Archivos con banderas: %d de %d\n', nnz(~cellfun(@isempty, T.banderas)), height(T));
