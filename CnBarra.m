function[w]=CnBarra(Cn, Cb, B, N)
% Función que calcula los costos reducidos
    w = Cn - Cb*(B^-1)*N;
end