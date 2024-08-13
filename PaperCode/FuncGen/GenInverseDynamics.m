function GenInverseDynamics(filename, E, P, varargin)


    p = inputParser;
    p.addOptional('mex', true);
    p.addOptional('return_Y', false);
    p.addOptional('help', true);
    p.addOptional('tol', 1e-9);

    p.parse(varargin{:});
    opt = p.Results;
    opt.alg = 'ID';
    
    DOF = size(P, 3);
    
    
    [vars, names, code] = GenSetupCode(E, P, opt);
    [vars, names, code] = CodeP(vars, names, code, opt);
    
    tauName = arrayfun(@(i) sprintf('tau[startInd+%d]', i), 0:DOF-1, 'UniformOutput', false);
    [tauCode, names, code] = GenTauCode(vars.E, vars.P, tauName, 'Y', opt.return_Y, names, code, opt);
    code.tau = tauCode;
        
    
    [funcDir, names.func, ext] = fileparts(filename);
    
    if opt.mex
        names.mexFunc = names.func;
        names.func = strcat(names.mexFunc, '_lib');
    end
    
    dir = fileparts(mfilename('fullpath'));
    
    allName = '';
    if opt.return_Y; allName = '_all'; end
    funcText = fileread( fullfile(dir, 'templates', sprintf('inverseDynamics%s.c', allName)));
    
    funcText = Replacements( funcText, vars, names, code, opt );
    
    functionPath = fullfile(funcDir, strcat(names.func, '.c'));
    fid = fopen(functionPath, 'wt');
    fprintf(fid, '%s', funcText);
    fclose(fid);
    
    headerText = fileread( fullfile(dir, 'templates', sprintf('inverseDynamics%s.h', allName)) );
    headerText = Replacements( headerText, vars, names, code, opt );

    headerPath = fullfile( funcDir, strcat(names.func, '.h') );
    fid = fopen( headerPath , 'wt');
    fprintf(fid, '%s', headerText);
    fclose(fid);
    
    
    if opt.mex
        
        mexText = fileread( fullfile(dir, 'templates', sprintf('inverseDynamics_mex%s.c', allName)) ); 
        mexText = Replacements( mexText, vars, names, code, opt );
        
        fid = fopen(fullfile( funcDir, strcat(names.mexFunc, '.c') ), 'wt');
        fprintf(fid, '%s', mexText);
        fclose(fid);
        
        fprintf("Compiling into mex file...\n");
        mex(filename, '-R2018a', '-outdir', funcDir, functionPath)
        fprintf("Done!\n");
        
        if opt.help
            helpText = fileread( fullfile(dir, 'templates', sprintf('inverseDynamics_help%s.m', allName)) );
            helpText = Replacements( helpText, vars, names, code, opt );
            help_filename = fullfile( funcDir, strcat(names.mexFunc, '.m') );
            fid = fopen(help_filename, 'wt');
            fprintf(fid, '%s', helpText);
            fclose(fid);
        end
        
    end
    
end