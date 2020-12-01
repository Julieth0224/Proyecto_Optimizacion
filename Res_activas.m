% Activas
function[Ak,I]=Res_activas(A,b,x)
    ax = A*x;
    I = [];
    Ak = [];
    [fil,col] = size(ax);
    for i = 1:fil
        if ax(i) == b(i)
            I = [I, i];
            fila_activa = A(i,:);
            Ak = [Ak;fila_activa];
        end
    end
end