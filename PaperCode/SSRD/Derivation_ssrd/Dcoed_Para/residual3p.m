function r = residual3p(A,x,b)
    %% Parse inputs
    assert(size(A, 1) == size(b, 1), "A and b do not match in size");
    assert(size(A, 2) == size(x, 1), "A and x do not match in size");
    assert(size(x, 2) == size(b, 2), "x and b do not match in size");

    %% Try to run mex
    c = PSDM.config();
    if coder.target('matlab') && c.use_mex
        try
            r = utilities.residual3p_mex(A, x, b);
            return;
        catch e
            disp(e);
        end
    end

   %% Actual code
   m = size(x, 2);
   N = size(A, 1);
   
   r = permute( ...
            utilities.iparfor( @(i) residual1D( A(i, :), x, b(i, :) ), ...
                           N, [1, m], false), ...
            [3, 2, 1]);
        
end

function r_i = residual1D(A_i, x, b_i)

    m = size(x, 2);
    r_i = zeros(1, m);

    for j = 1:m
        r_i(j) = utilities.dot3p(A_i, x(:, j), -b_i(j));
    end

end