function [tauCode, names, code] = makeTauCode( E, P, tauName, Yname, keepY, names, code, opt)
    
    DOF = size(P, 3);
    ell = size(P, 2);
    tol = 1e-11;
    P_mask = abs(P) > tol;
    m = size(E, 2);
    n = size(E, 1);
    
    if m == 0 || n == 0
        tauCode = sprintf('\t\t%s = 0.0;\n', tauName{1});
        return;
    end
    
  
    wasSkipped = zeros((ell*DOF), 1, 'logical');
    tau_j = cell(ell*DOF, 1);
    yCode = cell(ell*DOF, 1);
    
    c1 = 0; 
    c2 = 0; 
    
    for j = 1:DOF
        
        codeRowEls = cell(0,1);
        c4 = 0; 
        
        for k = 1:ell       
        
            c1 = c1+1;

            rowMask = P_mask(:, k, j);
            colMask = any(E(:, rowMask) > 0, 2);
            
            if ~any(rowMask, 1)
                wasSkipped(c1) = true;
                yCode{c1} = sprintf('%s[startIndY + %d] = 0.0;', Yname, (k-1)*DOF+j-1);
                continue;
            end
            
            c4 = c4 + 1;
            
            Su_t = E(colMask, rowMask);
            
            if isempty(Su_t)                
                coef = mat2str(P(rowMask, k, j));
                yCode{c1} = sprintf('%s[startIndY + %d] = %s;', Yname, (k-1)*DOF+j-1, coef);

                switch coef
                    case '1'
                        codeRowEls{c4} = sprintf('Theta[%d]', k-1);
                    case '-1'
                        codeRowEls{c4} = sprintf('-Theta[%d]', k-1);
                    otherwise
                        codeRowEls{c4} = sprintf('%s*Theta[%d]', coef, k-1);
                end
                
            else
                c2 = c2+1;

                names_t.gamma{1} = names.gamma{1}(colMask);
                names_t.gamma{2} = names.gamma{2}(colMask);

                names_t.y = cell(nnz(rowMask),1);
                c3 = 0;
                for i = 1:m
                    if rowMask(i)
                        c3 = c3+1;
                        names_t.y{c3} = mat2str(P(i, k, j));
                    end
                end

                subInd = sub2ind([DOF, ell], j, k);

                name = sprintf('%s_p%d', Yname, subInd);
                yCode{c1} = sprintf('%s[startIndY + %d] = %s;', Yname, (k-1)*DOF+j-1, name);

                [tau_j{c2}, names, code] = GenCode(Su_t, ...
                    name, ...
                    sprintf('%s_%d', Yname, c2),...
                    names_t, names, code, opt, 1);
           
                codeRowEls{c4} = sprintf('%s*Theta[%d]', name, k-1);
                
            end
            
        end
        
        if c4 > 0
            codeRows{j} = strrep( strjoin(codeRowEls, '+'), '+-', '-');
        else 
            codeRows{j} = '0.0';
        end
                
    end
        

    if c2 > 0
        tauCode = strjoin(regexprep(tau_j(1:c2), '\n+', '\n'), '\n\t\t');
    else
        tauCode = '';
    end
    
    for i = 1:DOF
        tauCode = sprintf('%s\n\t\t%s = %s;', tauCode, tauName{i}, codeRows{i});
    end
    
    if keepY
        tauCode = sprintf('%s\n\n\t\t%s', tauCode, strjoin(yCode, '\n\t\t'));
    end
       
end
