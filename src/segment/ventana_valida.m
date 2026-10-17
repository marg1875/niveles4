function ok = ventana_valida(malas, inicios, L)
%VENTANA_VALIDA  true si la ventana [inicio, inicio+L-1] no toca tramos malos.
%   ok = ventana_valida(malas, inicios, L), malas = k x 2 de marcar_artefactos
%   (vacío -> todas válidas). inicios puede ser un vector.
ok = true(size(inicios));
if isempty(malas), return, end
for i = 1:numel(inicios)
    a = inicios(i); b = a + L - 1;
    ok(i) = ~any(malas(:, 1) <= b & malas(:, 2) >= a);
end
end
