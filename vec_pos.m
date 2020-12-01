% Funcion que ve si todo el vector es psositivo
function[res, neg_pos]=vec_pos(x)
    res = true;
    neg_pos = 0;
    [fil,col] = size(x);
    neg = 0;
    j = 0;
    for i = 1:fil
        if x(i) <= 0
            j = j+1;
            if neg > x(i)
                neg = x(i);
                neg_pos = i;
            end
        end
    end
    if j > 0
        res = false;
    end
end