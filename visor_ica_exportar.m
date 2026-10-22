%% VISOR_ICA_EXPORTAR  Exporta señales y métricas ANTES y DESPUÉS del ICA.
%  Reproduce el preprocesamiento final (filtro -> ASR k=30 -> CAR -> z-score
%  -> blanqueo + ICA por sesión), igual que preprocesar_sesion.m, pero
%  guardando la señal justo antes del ICA, la señal limpia, las fuentes, la
%  matriz de mezcla A y la curtosis de cada componente.
%  Verifica que la señal limpia coincida con resultados_asr/senales.
%  -> resultados_asr/visor_ica/datos_ica.json (lo usa visor_ica.html)

raiz = fileparts(mfilename('fullpath'));
addpath(raiz, genpath(fullfile(raiz, 'src')));
setenv('NIVELES_RESULTADOS', 'resultados_asr');
cfg = config_niveles();
dirO = fullfile(cfg.rutas.resultados, 'visor_ica');
if ~isfolder(dirO), mkdir(dirO); end

fs = cfg.fs;
dec = 2;                                   % extractos a 125 Hz (señal <= 45 Hz)
dur_ext = 8;                               % segundos por extracto
bandas = [0.5 4; 4 8; 8 13; 13 30; 30 45];
esc = 20;                                  % extractos cuantizados a 0.05 z

T = listar_registros(cfg);
V = T(strcmp(T.estado, 'ok'), :);
[g, pac, mes, cond] = findgroups(V.paciente, V.mes, V.condicion);
F = disenar_filtros(cfg);
% VISOR_SESIONES="Px006_M1_FIS,Px007_M1_FIS" limita la exportación (vista rápida)
solo = strsplit(getenv('VISOR_SESIONES'), ',');
solo = solo(~cellfun(@isempty, solo));
if ~isempty(solo)
    dirO = fullfile(dirO, 'parcial');
    if ~isfolder(dirO), mkdir(dirO); end
end
P =coordenadas_canales(cfg.canales.nombres);

ses = {};
for s = 1:max(g)
    id = sprintf('%s_M%d_%s', pac{s}, mes(s), cond{s});
    if ~isempty(solo) && ~ismember(id, solo), continue; end
    filas = sortrows(V(g == s, :), 'nivel');
    rutas = filas.ruta.';
    nf = numel(rutas);
    % ---- igual que preprocesar_sesion (puentes = 'ninguno')
    Xs = cell(1, nf); fasr = zeros(1, nf);
    for j = 1:nf
        X = filtrar_eeg(cargar_registro(rutas{j}), F);
        E = asr_calibrar(X, cfg);
        [X, fasr(j)] = asr_procesar(X, E, cfg);
        X = referencia_promedio(X);
        Xs{j} = normalizar_zscore(X);
    end
    Xc = cell2mat(Xs);
    Xe = Xc(:, 1:cfg.ica.submuestreo:end);
    [~, Q, ~] = blanquear(Xe, cfg);
    Z = Q * (Xe - mean(Xe, 2));
    [W, ~, iinfo] = ica_gradiente_natural(Z, cfg);
    Wt = W * Q;
    Y = Wt * (Xc - mean(Xe, 2));
    [idx, curt] = limpiar_artefactos_ica(Y, cfg);
    lk = log(max(curt + 3, eps));
    zc = (lk - median(lk)) ./ max(1.4826 * median(abs(lk - median(lk))), eps);
    A = pinv(Wt);
    % varianza que aporta cada componente a los electrodos (z-score -> ~1 por canal)
    vy = var(Y, 0, 2);
    aporte = (A.^2) .* vy.';               % 16 x k

    % ---- verificación contra las señales guardadas
    Sg = load(fullfile(cfg.rutas.resultados, 'senales', [id '.mat']), 'X');
    err = 0;

    arch = cell(1, nf);
    for j = 1:nf
        pre = Xs{j};
        Yj = Wt * (pre - mean(Xe, 2));
        post = pre;
        if ~isempty(idx), post = pre - A(:, idx) * Yj(idx, :); end
        err = max(err, max(abs(double(Sg.X{j}(:)) - post(:))));
        a.archivo = nombre_corto(rutas{j});
        a.nivel = filas.nivel(j);
        a.frac_asr = fasr(j);
        a.n = size(pre, 2);
        a.pre = metricas(pre, fs, bandas);
        a.post = metricas(post, fs, bandas);
        % extracto: donde más pesa lo removido (o a mitad del registro)
        L = dur_ext * fs;
        if isempty(idx)
            c0 = round(a.n / 2);
        else
            r = sum((A(:, idx) * Yj(idx, :)).^2, 1);
            r = movmean(r, fs);
            [~, c0] = max(r(L/2+1 : end-L/2));
            c0 = c0 + L/2;
        end
        rg = (c0 - L/2 + 1) : (c0 + L/2);
        a.t0 = (rg(1) - 1) / fs;
        a.ext_pre = cuant(pre(:, rg(1:dec:end)), esc);
        if isempty(idx), a.ext_post = []; else, a.ext_post = cuant(post(:, rg(1:dec:end)), esc); end
        Ye = Yj(:, rg(1:dec:end));
        a.ext_fuentes = cuant(Ye ./ std(Yj, 0, 2), esc);
        arch{j} = a;
    end
    e.sesion = id;
    e.removidos = idx;
    e.curtosis = curt(:).';
    e.z_curtosis = zc(:).';
    e.iteraciones = iinfo.iteraciones;
    e.convergio = iinfo.convergio;
    e.A = A;
    e.aporte = aporte;
    e.err_vs_guardado = err;
    e.archivos = arch;
    ses{end+1} = e; %#ok<SAGROW>
    fprintf('%s: removidas=%s, iter=%d, |post - guardado| max = %.2e\n', id, mat2str(idx), iinfo.iteraciones, err);
end

D.canales = cfg.canales.nombres;
D.pos = P;
D.fs_ext = fs / dec;
D.escala = esc;
D.bandas = bandas;
D.nombres_bandas = {'delta','theta','alfa','beta','gamma'};
D.f_psd = metricas(zeros(16, 4 * fs), fs, bandas).f;
D.umbral_z = cfg.ica.curtosis_nsigma;
D.max_removidos = cfg.ica.max_removidos;
D.sesiones = ses;
fid = fopen(fullfile(dirO, 'datos_ica.json'), 'w', 'n', 'UTF-8');
fwrite(fid, jsonencode(D), 'char');
fclose(fid);
% visor autocontenido: la plantilla con los datos incrustados
html = fileread(fullfile(raiz, 'src', 'reportes', 'visor_ica_plantilla.html'));
html = strrep(html, '/*__DATOS__*/', fileread(fullfile(dirO, 'datos_ica.json')));
fid = fopen(fullfile(dirO, 'visor_ica.html'), 'w', 'n', 'UTF-8');
fwrite(fid, html, 'char');
fclose(fid);
fprintf('Listo: %s\n', fullfile(dirO, 'visor_ica.html'));

function m = metricas(X, fs, bandas)
% Métricas por canal sobre el registro completo.
[pxx, f] = pwelch(X.', 2 * fs, fs, 2 * fs, fs);
sel = f >= 0.5 & f <= 45;
m.f = f(sel).';
m.psd_db = round(10 * log10(pxx(sel, :).' + 1e-12), 2);   % 16 x nf
m.var = var(X, 0, 2).';
m.curtosis = (mean((X - mean(X, 2)).^4, 2) ./ var(X, 1, 2).^2 - 3).';
m.pico = max(abs(X), [], 2).';
bp = zeros(size(X, 1), size(bandas, 1));
for b = 1:size(bandas, 1)
    k = f >= bandas(b, 1) & f < bandas(b, 2);
    bp(:, b) = trapz(f(k), pxx(k, :)).';
end
m.banda = bp;                                   % potencia absoluta (z^2)
m.banda_rel = bp ./ sum(bp, 2);
C = corrcoef(X.');
m.corr = round(C, 3);
end

function q = cuant(X, esc)
q = int16(max(min(round(X * esc), 32767), -32767));
end

function n = nombre_corto(r)
[~, n] = fileparts(r);
end
