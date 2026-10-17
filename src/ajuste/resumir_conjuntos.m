function Rs = resumir_conjuntos(T)
%RESUMIR_CONJUNTOS  Exactitud media por conjunto (sesiones x clasificadores) y por clasificador.
[g, cj] = findgroups(T.conjunto);
Rs = table(cj, splitapply(@mean, T.acc, g), splitapply(@std, T.acc, g), splitapply(@mean, T.kappa, g), ...
    'VariableNames', {'conjunto','acc_media','acc_de','kappa_media'});
for c = unique(T.clasificador).'
    [g2, cj2] = findgroups(T.conjunto(strcmp(T.clasificador, c{1})));
    v = splitapply(@mean, T.acc(strcmp(T.clasificador, c{1})), g2);
    [~, ia, ib] = intersect(Rs.conjunto, cj2);
    Rs.(['acc_' c{1}]) = nan(height(Rs), 1);
    Rs.(['acc_' c{1}])(ia) = v(ib);
end
Rs = sortrows(Rs, 'acc_media', 'descend');
end
