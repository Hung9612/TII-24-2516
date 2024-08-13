function n2 = norm2(v, dim)
    
    if nargin < 2 || isempty(dim)
        dim = 1;
    end
    
    n2 = sum(v.^2, dim);
    
end