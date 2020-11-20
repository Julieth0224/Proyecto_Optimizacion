function[xa, Ib, In, N]=va_art(A, x, I)
% Creación de las variables artificiales necesarias para el problema
[fila, colum]=size(A);
s=length(x);
[fil, col]=size(I);
Ib=[];
In=[];
xa=[];
    for i=1:col
        l=0;
        for q=1:colum
            if A(:,q)==I(:,i)
                Ib = [Ib, q];
                l = l+1;
                break 
            end
        end
        if l == 0
            Ib = [Ib, s+1];
            xa = [xa, s+1];
            s = s+1;
        end
    end
    
    for d=x
        if not(ismember(d,Ib))
            In = [In, d];
        end
    end
N=[];
for i=1:length(In)
    N(:,i) = A(:,In(i));
        
end  
        
end
