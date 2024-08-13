function [vars, names, code] = GenSetupCode(E, P, opt)

    vars.E = E;
    vars.P = P;
    vars.DOF = size(P, 3);
    DOF = vars.DOF;
    vars.ell = size(P, 2);
    ell = vars.ell;
    
    if strcmp(opt.alg, 'ID')
        colMask = any(E > 0, 2);
    else
        colMask = any(E(1:(4*DOF), :) > 0, 2);
    end
    squareMask = any( E(colMask, :) > 1, 2);
    vars.colMask = colMask;
    vars.squareMask = squareMask;
    
    code.setup1 = '';
    code.setup2 = '';
    code.extra = '';
    code.tau = '';
    code.name_def = '';
    
    code.all = cell(0, 1);
    names.all = cell(0, 1);
    
    names.def = cell(0, 1);
    
    
    n = size(E, 1);
    ind = 1:n;
    ind1 = ind(colMask);
    ind2 = ind1(squareMask);
    
    names.gamma{1} = repelem({''}, n);
    names.gamma{2} = repelem({''}, n);
    
    for i = 1:nnz(colMask)
        names.gamma{1}{ ind1(i) } = sprintf('gi[%d]', ind1(i)-1);
    end
    
    code.gamma{2} = cell(nnz(squareMask), 1);
    for i = 1:nnz(squareMask)
        names.gamma{2}{ ind2(i) } = sprintf('gi2_p%d', i);
        code.gamma{2}{ i } = sprintf('gi[%d]*gi[%d]', ind2(i)-1, ind2(i)-1);
    end
    
    [code.gamma2, names] = assignVectors( 'gi2', code.gamma{2}, names, opt );
    code.setup2 = code.gamma2;

    if strcmp(opt.alg, 'FD')
        vars.accelMask_i = E((1:vars.DOF)+4*vars.DOF, :) > 0;
        vars.accelMask = any(vars.accelMask_i, 1);
        vars.Eind = E( :, ~ vars.accelMask);
    end
        
end