function[I,J]=itpl(L, xin ,h , Ng ,order)
% ITPL.M - Interpolation function for order 1-3
    I = zeros(Ng,1);
    
    in_meter = L*xin;
    l0 = floor(in_meter/h);
    al0 = in_meter/h - l0;

   
    if order == 1
        if l0-1 <= 0
            I(1)=1;
            J = 1/h*I;
        else
        I(l0) = 1;
        J = 1/h * I;
        end
    end

    if order == 2
        if l0-1 <= 0
            I(1)=1;
            J = 1/h*I;
        else
        I(l0) = (1-al0);
        I(l0 + 1) = al0;
        J = 1/h * I;
        end
    end

    if order == 3
        if l0-1 <= 0
            I(1)=1;
            J = 1/h*I;
        else
        I(l0 - 1) = (al0 * (al0 - 1) * (al0 - 2)) / -6;
        I(l0) = ((al0 - 1) * (al0 + 1) * (al0 - 2)) / 2;
        I(l0+1) = (al0 * (al0 + 1) * (al0 - 2)) / -2;
        I(l0 + 2) = (al0 * (al0 + 1) * (al0 - 1)) / 6;
        J = 1/h * I;
        end
    end

    I = I';
end