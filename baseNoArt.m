function[r]=baseNoArt(Ib, xa)
%Funcion que mira si existen variables artificiales en la base
r=0;
s=[];
    for i=1:length(Ib)
        for n=1:length(xa)
            if xa(n)==Ib(i)
                s = [s, n];
            end
        end   
    end
    if isempty(s)
        r=0;
    else 
        r=s(1);
    end

end
