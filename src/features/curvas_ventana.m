function f = curvas_ventana(V, cfg, baseH, lay)
%CURVAS_VENTANA  Curvas de Higuchi y DFA (y referencias) de una ventana.
%   f = curvas_ventana(V, cfg, baseH, lay), V (C x L). Fila con, en orden:
%     caracterizadores_limites(V, baseH)        familias Hurst de referencia
%     higuchi_dimension(V)  / dfa_exponente(V)  valores por defecto de Niveles3
%     higuchi_curva(V, lay.kmax)                L(k), k = 1..kmax
%     dfa_curva(V, lay.ns, orden, paso)         F(n) para cada orden x paso
%   (bloques de lay.bloques). Las curvas se guardan como reshape(C x K, 1, [])
%   (canal más rápido). De ellas salen todos los candidatos de la
%   calibración por regresión en un rango (caracteristicas_desde_curvas).
f = [caracterizadores_limites(V, baseH, cfg), ...
     higuchi_dimension(V, cfg).', dfa_exponente(V, cfg).', ...
     reshape(higuchi_curva(V, lay.kmax), 1, [])];
for a = lay.ordenes
    for b = 1:numel(lay.pasos)
        f = [f, reshape(dfa_curva(V, lay.ns, a, lay.pasos(b)), 1, [])]; %#ok<AGROW>
    end
end
end
