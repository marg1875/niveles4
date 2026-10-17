function P = coordenadas_canales(nombres)
%COORDENADAS_CANALES  Posiciones 3D aproximadas (esfera unitaria) del sistema 10-20/10-10.
%   P = coordenadas_canales(nombres) devuelve una matriz numel(nombres) x 3.
%   Coordenadas esféricas estándar (azimut, elevación en grados) de las
%   posiciones del montaje g.Nautilus del estudio. "Fp3" no es una posición
%   estándar; se asume Fp2 (simétrica de Fp1).
pos = struct( ...
  'Fp1', [108 -2],  'Fp3', [72 -2],  'Fp2', [72 -2], ...
  'F3',  [130 30],  'Fz',  [90 46],  'F4',  [50 30], ...
  'Cz',  [0 90],    'Pz',  [-90 46], 'P3',  [-130 30], 'P4', [-50 30], ...
  'P5',  [-140 13], 'P6',  [-40 13], 'P8', [-36 -2], ...
  'O1',  [-108 -2], 'Oz',  [-90 -2], 'T8', [0 -2],   'T7', [180 -2]);
P = zeros(numel(nombres), 3);
for i = 1:numel(nombres)
    ae = pos.(nombres{i});
    az = deg2rad(ae(1)); el = deg2rad(ae(2));
    P(i, :) = [cos(el) * cos(az), cos(el) * sin(az), sin(el)];
end
end
