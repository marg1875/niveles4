function [tau, RS] = rangors(v, particiones)

N = length(v);
i = 0;
m = N; 
tau = [];
RS = [];

while m >= (N / particiones) % dependiendo del numero de particiones a realizar
    R = [];
    S = [];
    j = 1;
    m = floor(N / (2 ^ i));
    while (j * m) <= N
        r = 0;
        M = 0;
        x = 0;
        V = 0; 
        r = v(1 + (j - 1) * m : j * m);
        M = mean(r);
        x = (r - M);               
        V = cumsum(x);
        R(j) = max(V) - min(V);
        S(j) = std(r);
        j = j + 1;
    end
    i = i + 1;
    tau(i) = m;
    RS(i) = mean(R / S);
end