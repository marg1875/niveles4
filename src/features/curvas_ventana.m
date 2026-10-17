function f = curvas_ventana(V, cfg, baseH, lay)
%CURVAS_VENTANA  Curvas de Higuchi y DFA (y referencias) de una ventana.
%   f = curvas_ventana(V, cfg, baseH, lay), V (C x L). Fila con, en orden:
%     caracterizadores_limites(V, baseH)        familias Hurst de referencia
%     higuchi_dimension(V)  / dfa_exponente(V)  valores por defecto de Niveles3
%     higuchi_curva(V, lay.kmax)                L(k), k = 1..kmax
%     dfa_curva(V, lay.ns, 1), dfa_curva(V, lay.ns, 2)   F(n) DFA1 y DFA2
%   Las curvas se guardan como reshape(C x K, 1, []) (canal más rápido).
%   De ellas salen todos los candidatos de la calibración por regresión en
%   un rango (caracteristicas_desde_curvas), sin recalcular nada.
f = [caracterizadores_limites(V, baseH, cfg), ...
     higuchi_dimension(V, cfg).', dfa_exponente(V, cfg).', ...
     reshape(higuchi_curva(V, lay.kmax), 1, []), ...
     reshape(dfa_curva(V, lay.ns, 1), 1, []), reshape(dfa_curva(V, lay.ns, 2), 1, [])];
end
