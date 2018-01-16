function [pcoef,x] = Ch1_Interp(f,a,b,n,method)
    if nargin<5
        method = 'uniform'
    end

    switch method
        case 'uniform',
            x=linspace(a,b,n);
            pcoef=polyfit(x,f(x),4);
    end
end