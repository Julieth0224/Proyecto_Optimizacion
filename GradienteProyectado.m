% Funcion que realiza el  algoritmo de gradiente proyectado
function[xmin]=GradienteProyectado(f, fgrad, A, epsilon, xk, b)
    xmin = [inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf; inf];
    [Ak, I] = Res_activas(A,b,xk);
    while norm(xmin-xk)> epsilon
        if xmin(1) ~= inf
            xk = xmin;
        end
        [fil, col] = size(Ak);
        [fili, coli] = size(I);
        P = eye(col)-Ak'*(Ak*Ak')^(-1)*Ak;
        gfxk = fgrad(xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11),xk(12),xk(13),xk(14),xk(15),xk(16),xk(17),xk(18),xk(19),xk(20),xk(21),xk(22),xk(23),xk(24),xk(25),xk(26),xk(27),xk(28));
        dk = -P*gfxk;
        if dk == zeros(size(dk))
            mu = -(Ak*Ak')^(-1)*Ak*gfxk;
            [res, neg] = vec_pos(mu);
            if res
                return;
            else
                temp = [];
                for i = 1:coli
                    if i ~= neg
                        temp = [temp, I(i)];
                    end
                end
                I = temp;
                [filii, colii] = size(I);
                Ak_new = [];
                for i = 1:colii
                    Ak_new = [Ak_new; A(I(i),:)];
                end
                Ak = Ak_new;
            end
        else
            syms alfa;
            falfa= xk + alfa*dk;
            v_res = A*falfa;
            [fiv,cov] = size(v_res);
            for q = 1:fiv
                v_res(q) = v_res(q)-b(q);
            end
            v_res_new = [];
            for q = 1:fiv
                if v_res(q) ~= 0
                    v_res_new = [v_res_new v_res(q)];
                end
            end
            [fivn,covn] = size(v_res_new);
            alfas = [];
            s = 0;
            for i = 1:covn
                s = solve(v_res_new(i),[alfa]);
                alfas = [alfas s];
            end
            if alfas >= 0
                alfa1 = min(alfas);
            elseif alfas < 0
                alfa1 = 0;
            else
                [filal, colal] = size(alfas);
                alfas_pos = [];
                for q = 1:colal
                    if alfas(q)>=0
                        alfas_pos = [alfas_pos alfas(q)];
                    end
                end
                alfa1 = min(alfas_pos);
            end
            vec = xk + alfa*dk;
            fun = f(vec(1), vec(2),vec(3), vec(4),vec(5), vec(6),vec(7), vec(8),vec(9), vec(10),vec(11), vec(12),vec(13), vec(14),vec(15), vec(16),vec(17), vec(18),vec(19), vec(20),vec(21), vec(22),vec(23), vec(24),vec(25), vec(26),vec(27), vec(28));
            fun2 =@(w) subs(fun, alfa, w);
            [alfa2, t] = MetodoAureo(fun2,0,alfa1,false);
            xmin = xk + alfa2*dk;
        end
    end  
end