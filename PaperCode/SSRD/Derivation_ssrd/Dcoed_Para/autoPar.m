function doPar = autoPar

    c = config;
    if coder.target('matlab')
        doPar = c.use_par && ~isempty(gcp('nocreate'));
    else
        doPar = c.use_par;
    end
    
end
