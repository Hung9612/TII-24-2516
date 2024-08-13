function X = pinv(A, varargin)

    X = pinv(A, varargin{:});
    
    r = utilities.residual3p(A, X, eye(size(A, 1)));  
    
    X = X - X*r;

end