function vprint(verbose, varargin)
    
    if verbose
        fprintf(varargin{:});
    end
end