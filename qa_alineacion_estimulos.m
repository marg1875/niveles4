%% QA_ALINEACION_ESTIMULOS  Evidencia de la temporización uniforme de stims.
%
%   Los .mat no traen marcas temporales: stims es la secuencia de etiquetas
%   [0 0 0 0 L L L L] x 10. El pipeline asume que las 80 entradas se
%   reparten uniformemente sobre EEG_all. Este script documenta el supuesto:
%   para cada archivo pliega la potencia de artefactos sobre el ciclo de 8
%   entradas (4 de reposo + 4 de tarea):
%     - LF (0.5-3 Hz): movimiento / parpadeo
%     - HF (30-100 Hz): actividad muscular (EMG contaminante)
%   En la ejecución física (FIS) se espera más potencia de artefactos en la
%   mitad de tarea. contraste = media(tarea) - media(reposo) del log-perfil.
%   No hay canal de EMG en los archivos: los artefactos EEG son el proxy.
%   -> resultados/qa/alineacion.csv, resultados/qa/alineacion.png

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
cfg = config_niveles();
dirQ = fullfile(cfg.rutas.resultados, 'qa');
if ~isfolder(dirQ), mkdir(dirQ); end

T = listar_registros(cfg);
T = T(strcmp(T.estado, 'ok'), :);
fs = cfg.fs;
[bh, ah] = butter(4, [30 100]/(fs/2));
[bl, al] = butter(2, [0.5 3]/(fs/2));
n = height(T);
perfLF = zeros(n, 8); perfHF = zeros(n, 8);
for k = 1:n
    [X, st] = cargar_registro(T.ruta{k});
    X = X - mean(X, 2);
    N = size(X, 2); K = numel(st);
    bordes = round((0:K) * N / K);
    Ph = mean(filtfilt(bh, ah, X.').^2, 2).';
    Pl = mean(filtfilt(bl, al, X.').^2, 2).';
    ph = zeros(1, K); pl = zeros(1, K);
    for j = 1:K
        ph(j) = median(Ph(bordes(j)+1:bordes(j+1)));
        pl(j) = median(Pl(bordes(j)+1:bordes(j+1)));
    end
    fh = mean(reshape(log(ph), 8, []), 2).'; perfHF(k,:) = fh - mean(fh);
    fl = mean(reshape(log(pl), 8, []), 2).'; perfLF(k,:) = fl - mean(fl);
end
contraste_LF = mean(perfLF(:,5:8), 2) - mean(perfLF(:,1:4), 2);
contraste_HF = mean(perfHF(:,5:8), 2) - mean(perfHF(:,1:4), 2);
Q = [T(:, {'archivo','paciente','mes','condicion','nivel'}), ...
     array2table([perfLF, perfHF], 'VariableNames', [compose('LF_%d', 1:8), compose('HF_%d', 1:8)]), ...
     table(contraste_LF, contraste_HF)];
writetable(Q, fullfile(dirQ, 'alineacion.csv'));

f = figure('Visible', 'off', 'Position', [100 100 1100 800]);
subplot(1,2,1); imagesc(perfLF); colorbar; title('LF 0.5-3 Hz (log, centrado)');
xlabel('Entrada del ciclo (1-4 reposo | 5-8 tarea)'); yticks(1:n); yticklabels(T.archivo);
set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 6); xline(4.5, 'w', 'LineWidth', 2);
subplot(1,2,2); imagesc(perfHF); colorbar; title('HF 30-100 Hz (log, centrado)');
xlabel('Entrada del ciclo (1-4 reposo | 5-8 tarea)'); yticks(1:n); yticklabels(T.archivo);
set(gca, 'TickLabelInterpreter', 'none', 'FontSize', 6); xline(4.5, 'w', 'LineWidth', 2);
exportgraphics(f, fullfile(dirQ, 'alineacion.png'), 'Resolution', 150);
close(f);
fis = strcmp(T.condicion, 'FIS');
fprintf('Contraste tarea-reposo (FIS): LF %+.3f (positivo en %d/%d), HF %+.3f (positivo en %d/%d)\n', ...
    mean(contraste_LF(fis)), nnz(contraste_LF(fis) > 0), nnz(fis), ...
    mean(contraste_HF(fis)), nnz(contraste_HF(fis) > 0), nnz(fis));
fprintf('Resultados en %s\n', dirQ);
