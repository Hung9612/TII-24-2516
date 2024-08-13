function B = vertStack(A, fromDimension, squeezeResult)

    if nargin < 2 || isempty(fromDimension)
        fromDimension = 3;
    end
    squeezeResult = (nargin > 2) && squeezeResult;
    
    assert(fromDimension ~= 1, "Can't stack from dimension 1!");
    
    [n, m, p] = size(A);
        
    if fromDimension == 3
        B = reshape( permute(A, [1 3 2]), [n*p, m]);
    elseif fromDimension == 2
        if squeezeResult
            B = reshape( A, [n * m, p]);
        else
            B = reshape( A, [n * m, 1, p]);
        end
    else
        error("Invalid fromDimension!");
    end
    
end