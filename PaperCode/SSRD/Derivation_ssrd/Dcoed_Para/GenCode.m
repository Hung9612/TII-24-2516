function [setCode, names, code] = GenCode(S, name, name_prefix, names_t, names, code, opt, depth)
    

    n = size(S, 1);
    m = size(S, 2);
    
    if m < 1
        setCode = sprintf('\t\t%s = 0;\n', name);
        names = PSDM.fgen.addToNames(name, names);
        return;
    end
    if n < 1
        Psum = sum(cellfun(@str2num, names_t.y));
        setCode = sprintf('%s = %s;\n', name, mat2str(Psum));
        names = add2Names(name, names);
        return;
    end

    s = 1;
    
    S1 = S(1:s, :);
    S2 = S(s+1:end, :);

    Su1 = unique(S1', 'rows', 'stable')';
    Su2 = unique(S2', 'rows', 'stable')';
    
    m1 = size(Su1, 2);
    m2 = size(Su2, 2);
    
    if n == 1 || all(all( Su2 == 0, 1), 2)
        name_i = name;
        final = true;
    else
        name_i = sprintf('%s_d%d', name_prefix, depth);
        final = false;
    end
    
    c = 0;
    yCode = cell(0,0);
    yNames = cell(m2, 1);
    codeNameInd = zeros(m2, 1);
    for i = 1:m2
        
        s2 = Su2(:, i);
        
        s1_ind = find( all( S2 == s2, 1) );
        m1i = numel(s1_ind);
        
        s1 = S1( :, s1_ind );

        s1_empty = all(s1 == 0, 1);
        s1_single = m1i == 1;
        
        if ~s1_single || any(~s1_empty) || final
            
            forwardFlag = true;
            
            terms = cell(m1i,1);
            
            for k = 1:m1i
                if s1_empty(k)
                    terms{k} = sprintf('%s', names_t.y{s1_ind(k)});
                else
                    switch names_t.y{s1_ind(k)}
                        case '1'
                            terms{k} = sprintf('%s', names_t.gamma{ s1(1,k) }{ 1 });
                        case '-1'
                            forwardFlag = false;
                            terms{k} = sprintf('-%s', names_t.gamma{ s1(1,k) }{ 1 });
                        otherwise
                            forwardFlag = false;
                            terms{k} = sprintf('%s*%s', ...
                                names_t.y{s1_ind(k)}, names_t.gamma{ s1(1,k) }{ 1 });
                    end
                end
            end
            
            if m1i == 1 && forwardFlag && ~final
                yNames{i} = terms{1};
            else
                code_i = strrep( strjoin(terms, '+'), '+-', '-');

                name_already = findEquivNames(code_i, names, code);
                
                if isempty(name_already)
                    c = c+1;
                    yCode{c} = code_i;
                    codeNameInd(c) = i;
                    if final
                        yNames{i} = name_i;
                    else
                        yNames{i} = sprintf('%s_p%d', name_i, c);
                    end
                    
                    [code, names] = addNameCode(code_i, yNames{i}, code, names);
                    
                else
                    
                    yNames{i} = name_already;
                    if final
                        c = c+1;
                        yCode{c} = yNames{i};
                    end
                    
                end
            end
            
        else
            yNames{i} = names_t.y{s1_ind(1)};
            
        end
                
    end
    
    if c>0
        if final
            setCode = sprintf('%s = %s;', name_i, yCode{1} );
            names = add2Names(name_i, names);
        else
            [setCode, names] = assignVectors( name_i, yCode, names, opt );
        end
    else
        setCode = '';
    end
    
    if ~ final
        names_s = names_t;
        names_s.y = yNames;

        names_s.gamma{1} = names_t.gamma{1}(2:end);
        names_s.gamma{2} = names_t.gamma{2}(2:end);
        
        [setCode_s, names, code] = GenCode(Su2, name, name_prefix, names_s, names, code, opt, depth+1);
        
        setCode = sprintf('%s\n\n\t\t%s', setCode, setCode_s);
        
    end
    
end


function [code, names] = addNameCode(codeAdd, nameAdd, code, names)
    N = size(code.all, 1);
    code.all{N+1, 1} = codeAdd;
    names.all{N+1, 1} = nameAdd;
end