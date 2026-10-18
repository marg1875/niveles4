function f = caract_curvas(V, cfg, lay)
%CARACT_CURVAS  Curvas de escalamiento y valores de todas las familias de una ventana.
%   f = caract_curvas(V, cfg, lay), V (C x L) -> fila de lay.ancho valores,
%   en el orden de lay.bloques (caract_registro). Cada bloque es una matriz
%   C x K guardada como reshape(., 1, []) (canal más rápido). Los candidatos
%   (pendientes en rangos, valores con parámetros) salen después con
%   caract_desde_curvas, sin recalcular nada.
p = cfg.complejidad;
f = zeros(1, lay.ancho);
mf = struct();                                     % MFDFA: un cálculo por orden
for b = lay.bloques
    a = b.arg;
    switch b.nombre
        case 'ref'
            M = [reshape(caracterizadores_limites(V, a.baseH, cfg), [], numel(a.baseH)), ...
                 higuchi_dimension(V, cfg), dfa_exponente(V, cfg)];
        case 'hig', M = higuchi_curva(V, b.K);
        case 'dma', M = dma_curva(V, a.ns);
        case 'agv', [M, ~] = agregados_curva(V, a.ms);
        case 'abm', [~, M] = agregados_curva(V, a.ms);
        case 'psd', M = psd_curva(V, cfg.fs, a.nfft);
        case 'box', M = boxcount_curva(V, a.ss);
        case 'dir', M = dimensiones_directas(V);
        case 'samp', M = sampen_rejilla(V, p.se_m, p.se_r, p.se_tau);
        case 'fuz', M = fuzzyen_rejilla(V, 2, p.fe_r, p.fe_tau);
        case 'pe', M = permen_rejilla(V, p.pe_d, p.pe_tau);
        case 'lz', M = lzc_rejilla(V, p.lz_tau);
        case 'mse', M = mse_ci(V, p.mse_S, 2, 0.15);
        otherwise
            if startsWith(b.nombre, 'dfa')
                M = dfa_curva(V, a.ns, a.orden, a.paso);
            elseif startsWith(b.nombre, 'mf')
                k = sprintf('o%d', a.orden);
                if ~isfield(mf, k), mf.(k) = mfdfa_curva(V, a.ns, a.orden, p.mf_q); end
                M = mf.(k)(:, :, a.iq);
            elseif startsWith(b.nombre, 'wav')
                M = wavelet_curva(V, a.wname, a.J);
            end
    end
    f(b.o + (1:b.K * lay.C)) = reshape(M, 1, []);
end
end
