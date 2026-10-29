function P = coordenadas_canales(electrodos)
%COORDENADAS_CANALES  Posiciones 3D aproximadas (esfera unitaria) del sistema 10-20/10-10.
%   P = coordenadas_canales(electrodos) devuelve una matriz numel(electrodos) x 3.
%   Coordenadas esféricas estándar (azimut, elevación en grados) de las
%   posiciones del montaje g.Nautilus PRO de 16 canales. Recibe los nombres
%   reales (cfg.canales.electrodos), no los identificadores históricos.
pos = struct( ...
  'Fp1', [108 -2],  'Fp2', [72 -2], ...
  'F3',  [130 30],  'Fz',  [90 46],  'F4',  [50 30], ...
  'T7',  [180 -2],  'C3',  [180 46], 'Cz',  [0 90],   'C4',  [0 46],  'T8', [0 -2], ...
  'P3',  [-130 30], 'Pz',  [-90 46], 'P4',  [-50 30], ...
  'PO7', [-126 -2], 'PO8', [-54 -2], 'Oz',  [-90 -2]);
P = zeros(numel(electrodos), 3);
for i = 1:numel(electrodos)
    ae = pos.(electrodos{i});
    az = deg2rad(ae(1)); el = deg2rad(ae(2));
    P(i, :) = [cos(el) * cos(az), cos(el) * sin(az), sin(el)];
end
end
