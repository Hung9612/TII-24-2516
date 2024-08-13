function val = rateRound(in, tol, isSym)

    if nargin < 3
        isSym = false;
    end
    
    if nargin < 2
        tol = 1e-8;
    end
    
    [N, D] = rat(in, tol);

    if isSym
        val = sym(N)./D;
    else
        val = N./D;
    end

end
