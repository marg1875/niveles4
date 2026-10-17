function q = percentil95(x)
%PERCENTIL95  Percentil 95 (empírico) sin Statistics Toolbox.
x = sort(x(:)); q = x(max(1, ceil(0.95 * numel(x))));
end
