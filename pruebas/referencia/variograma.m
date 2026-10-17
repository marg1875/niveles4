function [vario, alt] = variograma(v, tam)

var = 0;
alt(1)=0;
% vario(1)=0;
for i = 2 : 1 : (tam/2)     %empieza en 2 porque se mueven los valores
    k = 1;               %puede ser modificado para ver todo el variograma
    vj = [];            %hay controlar este porque cuenta la cantidad de sumas que se hacen
    vk = [];
    var = 0;
    for j = 1 + i : 1 : tam
        vj(k) = v(j);
        vk(k) = v(k);
        var = var + (v(j) - v(k))^ 2;
%         var = var + (abs(v(j) - v(k)))^(1/2);
        k = k + 1;
    end
    if k > (24/100)*tam
        vario(i) = (var / (2 * (k-1)));% - (mean(vj) - mean(vk)) ^ 2;
%         vario(i) = ((0.457+0.494/(k-1)+0.494/((k-1)^2))^(-1)*((var*(1/(k-1)))^4))/2;
        %vario(i) = var / (2 * (tam - i)) - (mean(vj) - mean(vk)) ^ 2;
        %error de formula
        alt(i) = i - 1;
    else
        break
    end
end