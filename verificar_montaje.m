%VERIFICAR_MONTAJE  Contrasta dos asignaciones de electrodos a las filas de
%   EEG_all con las propias señales (2026-10-28):
%   - la lista histórica (cfg.canales.nombres) y
%   - el montaje estándar del g.Nautilus PRO (cfg.canales.electrodos).
%   Pruebas, con todos los registros que traen EEG_all:
%   1) Parpadeo: pendiente de cada canal sobre la media de las filas 1-2 en
%      0.5-4 Hz; debe decrecer de adelante hacia atrás.
%   2) Con referencia promedio (4-30 Hz), la correlación entre canales debe
%      decrecer con la distancia entre electrodos.
%   Uso: matlab -batch "verificar_montaje"
cfg = config_niveles();
addpath(genpath(fullfile(cfg.rutas.raiz, 'src')));
fs = cfg.fs;
[b, a] = butter(4, [4 30] / (fs / 2));
[bd, ad] = butter(4, [0.5 4] / (fs / 2));
L = dir(fullfile(cfg.rutas.datos, '**', '*.mat'));
Z = []; B = [];
for k = 1:numel(L)
    S = load(fullfile(L(k).folder, L(k).name));
    if ~isfield(S, 'EEG_all'), continue, end
    X = double(S.EEG_all);
    if size(X, 1) ~= 16, X = X.'; end
    X = X(:, 2 * fs:end).';                          % sin los primeros 2 s
    Y = filtfilt(b, a, X); Y = Y - mean(Y, 2);
    Z = cat(3, Z, atanh(min(corr(Y), 0.999)));
    D = filtfilt(bd, ad, X); x = mean(D(:, 1:2), 2);
    B(end + 1, :) = (D.' * x / (x.' * x)).'; %#ok<AGROW>
end
R = tanh(median(Z, 3)); iu = find(triu(true(16), 1)); bm = median(B);
fprintf('Registros con EEG_all: %d\n', size(B, 1));
% Ubicación de la lista histórica: posiciones 10-10 que no están en el g.Nautilus PRO
extra = struct('Fp3', [72 -2], 'P5', [-140 13], 'P6', [-40 13], 'P8', [-36 -2], 'O1', [-108 -2]);
listas = {cfg.canales.nombres, cfg.canales.electrodos};
titulos = {'Lista histórica', 'g.Nautilus PRO'};
for m = 1:2
    nom = listas{m};
    XY = zeros(16, 3);
    for c = 1:16
        if isfield(extra, nom{c}), ae = deg2rad(extra.(nom{c}));
            XY(c, :) = [cos(ae(2)) * cos(ae(1)), cos(ae(2)) * sin(ae(1)), sin(ae(2))];
        else, XY(c, :) = coordenadas_canales(nom(c)); end
    end
    Dm = squareform(pdist(XY));
    fprintf('\n== %s\n', titulos{m});
    fprintf('  Parpadeo: Spearman(pendiente, anterioridad), filas 3-16 = %.3f\n', ...
        corr(bm(3:16).', XY(3:16, 2), 'type', 'Spearman'));
    fprintf('  Referencia promedio: Spearman(correlación, cercanía) = %.3f\n', ...
        corr(R(iu), -Dm(iu), 'type', 'Spearman'));
end
fprintf('\nPendiente del parpadeo por fila (de mayor a menor):\n');
[~, ix] = sort(bm, 'descend');
for c = ix
    fprintf('  fila %2d  histórico %-4s  g.Nautilus %-4s  %.2f\n', c, cfg.canales.nombres{c}, cfg.canales.electrodos{c}, bm(c));
end
