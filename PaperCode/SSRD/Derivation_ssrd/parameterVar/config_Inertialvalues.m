function [numX,theta] = config_Inertialvalues(n,symtheta)
%CONFIG_INERTIALVALUE 
%
%
%   Output:
%   - X:    A DOF x 10 matrix of the mass properties of the
%           manipulator:
%           [Ixx_1  Iyy_1  Izz_1 Ixy_1 Ixz_1 Iyz_1 rcx_1 rcy_1 rcz_1 m_1;
%             :      :      :      :    :     :      :      :     :    :
%           [Ixx_n  Iyy_n  Izz_n Ixy_n Ixz_n Iyz_n rcx_n rcy_n rcz_n m_n];
%           If you include a DOF+1 row to this, it will assume that
%           this is the tool effect and combine it into the
%           equations.

if n==2
    X = [2   0.125 0 0 0 0 0 0 0 0.002;...
         2.5 0.125 0 0 0 0 0 0 0 0.002];
    numX=[X(:,1).',reshape(X(:,2:4).',[3*n,1]).',reshape(X(:,5:end).',[6*n,1]).'];
    symX=cal_Yd(n);
    theta=eval(subs(symtheta,num2cell(symX),num2cell(numX)));
    
elseif n==4
    X = [6.25  209.6  0      44.7   0.0104  0.1355 0.1413  0.0000  0.0056 -0.0000;...
        9.49  -166.2  0     -146.4  0.0494  0.1553 0.1336 -0.0000 -0.0129  0.0001;...
        0.40    0     0     -175.1  0.0069  0.0069 0.0000  0.0000  0.0000 -0.0000;...
        1.10   4.1   -12.4  -16.8   0.0023  0.0004 0.0025 -0.0000  0.0000  0.0001];
    numX=[X(:,1).',reshape(X(:,2:4).',[3*n,1]).',reshape(X(:,5:end).',[6*n,1]).'];
    symX=cal_Yd(n);
    theta=eval(subs(symtheta,num2cell(symX),num2cell(numX)));
elseif n==6
    X = [7.60  -0.0376  0.0435 -0.0076  0.0430 0.0449 0.0490  0.0125  0.0022 -0.0024;...
         2.80   0.1554 -0.0062  0.1148  0.0060 0.0393 0.0415 -0.0004 -0.0012  0.0001;...
         7.00   0.0450 -0.0070  0.0040  0.0278 0.0375 0.0309  0.0010  0.0005 -0.0002;...
         2.90  -0.0010  0.1021 -0.0028  0.0138 0.0058 0.0123 -0.0010  0.0002  0.0002;...
         0.94   0.0001 -0.0300  0.0106  0.0009 0.0010 0.0006  0.0001  0.0010  0.0001;...
         3.00  -0.0020 -0.0190  0.0670  0.0045 0.0050 0.0060  0.0001 -0.0003  0.0004];
    numX=[X(:,1).',reshape(X(:,2:4).',[3*n,1]).',reshape(X(:,5:end).',[6*n,1]).'];

    symX=cal_Yd(n);
    theta=eval(subs(symtheta,num2cell(symX),num2cell(numX)));
    
elseif n==7
    X = [3.95 -0.35 1.600 -31.40 0.0046  0.0045 0.0003  0.0000  0.0000 -0.0000;...
         4.50 -7.7  166.8 -3.600 0.0003  0.0001 0.0004 -0.0000  0.0000  0.0000;...
         2.45 -2.2 -34.90 -26.50 0.0022  0.0022 0.0007  0.0001  0.0001  0.0001;...
         2.61  0.2 -52.70  38.20 0.0384  0.0114 0.0499 -0.0009 -0.0011  0.0011;...
         3.41  0.1 -2.400 -211.3 0.0028  0.0028 0.0001  0.0000  0.0000  0.0000;...
         3.39  0.5  20.20  27.50 0.0005  0.0028 0.0023 -0.0001  0.0000  0.0000;...
         3.00 -2.0 -19.00  67.00 0.0045  0.0050 0.0060  0.0001 -0.0003  0.0004];
    numX=[X(:,1).',reshape(X(:,2:4).',[3*n,1]).',reshape(X(:,5:end).',[6*n,1]).'];
    symX=cal_Yd(n);
    theta=eval(subs(symtheta,num2cell(symX),num2cell(numX)));
    
else
    sprintf('This kind of robot should be config in function config_DH.m')
end

end

