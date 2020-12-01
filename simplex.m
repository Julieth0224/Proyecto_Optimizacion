function[xf]=simplex(A, c, b, max, mode)
% Funcion que realiza el metodo simplex
if max
    %Si se pide el maximo 
    c = -c;
end

x=[];
[fi,co]=size(A);
for r=1:co
    x = [x, r];
end
I=eye(fi);
B=I;
[xa, Ib, In, N]=va_art(A, x, I);
xb=b;
xf = [];

z=baseNoArt(Ib, xa);

% Modo silent 
if mode == false
    % Realizamos el metoso de dos fases: 
    % Fase 1
    while z ~= 0
        if det(B)==0
            disp('Matriz no invertible.')
            return;
        end
        [f,col]=size(N);
        cba=[];
        for w=Ib
            if ismember(w,xa)
                cba = [cba, 1];
            else
                cba = [cba, 0];
            end
        end
        cost = CnBarra(zeros(1,col), cba, B, N);
        if min(cost)>=0
            disp('No tiene solución básica factible.')
            return;
        end
    
        for i=1:length(cost)
            if min(cost)==cost(i)
                pos_entra = i;
            end
        end
        entra=In(pos_entra);
        Y = Yk(B, A(:,entra));
        if all(Y<0)
            disp('No tiene optimo finito.')
            return;
        end
        [e, ei]=razon(Y, xb);
        for i=1:length(e)
            if min(e)==e(i)
                pos_sale = ei(i);
                break
            end
        end
        sale=Ib(pos_sale);
        for i=1:length(Ib)
            if Ib(i)==sale
                Ib(i)=entra;
                break
            end
        end 
    
        B(:,pos_sale)=N(:,pos_entra);
        for i=1:length(xa)
            if sale == xa(i)
                In(pos_entra)=[];
                N(:,pos_entra)=[];
            else
                In(pos_entra)=sale;
                N(:,pos_entra)=A(:,sale);
            end
        end
        z=baseNoArt(Ib, xa);
        xb=(B^-1)*(b);
    end
    
    % Fase 2
    cb=[];
    cn=[];
    cb=c(Ib);
    cn=c(In);
    if det(B)==0
        disp('Matriz no invertible.')
        return;
    end
    cost2= CnBarra(cn, cb, B, N);

    while min(cost2) < 0
        if det(B)==0
            disp('Matriz no invertible.')
            return;
        end
        for i=1:length(cost2)
            if min(cost2)==cost2(i)
                pos_entra2 = i;
            end
        end
        entra2=In(pos_entra2);
        Y2= Yk(B, A(:,entra2));

        if all(Y2<0)
            disp('No tiene optimo finito.')
            return;
        end
        [e2, ei2]=razon(Y2, xb);
        for i=1:length(e2)
            if min(e2)==e2(i)
                pos_sale2 = ei2(i);
                break;
            end
        end
        sale2=Ib(pos_sale2);
        for i=1:length(Ib)
            if Ib(i)==sale2
                Ib(i)=entra2;
                break
            end
        end 
        In(pos_entra2)=sale2;
        B=[];
        for t=Ib
            B=[B, A(:,t)];
        end
        N=[];
        for f=In
            N=[N, A(:,f)];
        end
        cb=c(Ib);
        cn=c(In);
        cost2= CnBarra(cn, cb, B, N);
        xb=(B^-1)*(b);
        if xb == zeros(length(xb),1)
            disp('Infinitas soluciones factibles.')
            return;
        end
    end

    for s=x
        if ismember(s,Ib)
            xf = [xf, xb(Ib==s)];
        else
            xf = [xf, 0];
        end
    end
    if xf == zeros(1, length(xf))
        disp('Infinitas soluciones factibles')
        return;
    end


% Modo verbose 
elseif mode == true
    % Realizamos el metoso de dos fases: 
    % Fase 1
    
    disp('Indice de las variables básicas: ')
    disp(Ib)
    
    while z ~= 0
        if det(B)==0
            disp('Matriz no invertible.')
            return;
        end
        [f,col]=size(N);
        cba=[];
        for w=Ib
            if ismember(w,xa)
                cba = [cba, 1];
            else
                cba = [cba, 0];
            end
        end
        cost = CnBarra(zeros(1,col), cba, B, N);
        if min(cost)>=0
            disp('No tiene solución básica factible.')
            return;
        end
    
        for i=1:length(cost)
            if min(cost)==cost(i)
                pos_entra = i;
            end
        end
        entra=In(pos_entra);
        Y = Yk(B, A(:,entra));
        if all(Y<0)
            disp('No tiene optimo finito.')
            return;
        end
        [e, ei]=razon(Y, xb);
        for i=1:length(e)
            if min(e)==e(i)
                pos_sale = ei(i);
                break
            end
        end
        sale=Ib(pos_sale);
        for i=1:length(Ib)
            if Ib(i)==sale
                Ib(i)=entra;
                break
            end
        end 
    
        B(:,pos_sale)=N(:,pos_entra);
        for i=1:length(xa)
            if sale == xa(i)
                In(pos_entra)=[];
                N(:,pos_entra)=[];
            else
                In(pos_entra)=sale;
                N(:,pos_entra)=A(:,sale);
            end
        end
        z=baseNoArt(Ib, xa);
        xb=(B^-1)*(b);
        
        disp('Variable que entra: ')
        disp(entra)
        disp('Variable que sale: ')
        disp(sale)
        
        disp('Indice de las variables básicas: ')
        disp(Ib)
        
    end
    
    % Fase 2
    cb=[];
    cn=[];
    cb=c(Ib);
    cn=c(In);
    if det(B)==0
        disp('Matriz no invertible.')
        return;
    end
    cost2= CnBarra(cn, cb, B, N);

    while min(cost2) < 0
        if det(B)==0
            disp('Matriz no invertible.')
            return;
        end
        for i=1:length(cost2)
            if min(cost2)==cost2(i)
                pos_entra2 = i;
            end
        end
        entra2=In(pos_entra2);
        Y2= Yk(B, A(:,entra2));

        if all(Y2<0)
            disp('No tiene optimo finito.')
            return;
        end
        [e2, ei2]=razon(Y2, xb);
        for i=1:length(e2)
            if min(e2)==e2(i)
                pos_sale2 = ei2(i);
                break;
            end
        end
        sale2=Ib(pos_sale2);
        for i=1:length(Ib)
            if Ib(i)==sale2
                Ib(i)=entra2;
                break
            end
        end 
        In(pos_entra2)=sale2;
        B=[];
        for t=Ib
            B=[B, A(:,t)];
        end
        N=[];
        for f=In
            N=[N, A(:,f)];
        end
        cb=c(Ib);
        cn=c(In);
        cost2= CnBarra(cn, cb, B, N);
        xb=(B^-1)*(b);
        if xb == zeros(length(xb),1)
            disp('Infinitas soluciones factibles.')
            return;
        end
        tam=(cb)*(xb);
        
        disp('Variable que entra: ')
        disp(entra2)
        disp('Variable que sale: ')
        disp(sale2)
        
        disp('Indice de las variables básicas: ')
        disp(Ib)
        
        disp('Tamaño del paso: ')
        disp(tam)
    end

    for s=x
        if ismember(s,Ib)
            xf = [xf, xb(Ib==s)];
        else
            xf = [xf, 0];
        end
    end
    
    if xf == zeros(1, length(xf))
        disp('Infinitas soluciones factibles')
        return;
    end
    
else
    disp('Modo no aceptado.')
    return;
end



