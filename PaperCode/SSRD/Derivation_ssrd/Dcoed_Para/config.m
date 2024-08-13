function c = config()
    c.use_mex = true;    
    c.do_reprojection_tests = coder.target('matlab');
    
    c.use_par = false;
    
end
