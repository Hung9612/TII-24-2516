function [DH,isRevo] = config_DH(n)
%CONFIG_DH configurate the robot DH parameter table
% Input:
%   - n: DoF of robot. You can change the value of their DH values of the
%   specified DoF to meet your requirement.
%
% Outputs:
%   - DH_ext: a DOFx4-6 array of the DH params in the following order: 
%             [a_1  alpha_1    d_1   theta_1    lt_1    q_sign_1;
%              :       :        :       :         :         :     
%              an   alpha_n    d_n   theta_n    lt_n    q_sign_n];
%           lt stands for "Link types", and should be true if joint is a
%           prismatic joint, false otherwise. If column 5 is missing, will
%           assume revolute (0).
%           q_sign is a link direction number (-1 or 1). Allows you to
%           change the sign of the joint variable. If column 6 is
%           missing, will assume 1.
if n==2
    a1 = 250e-3;
    a2 = 0;
        % 
    DH = [a1,       0,        0,      0,      0,        1;
          a2,       pi,       0,      0,      0,        1];
           % 
    isRevo=[0,0,0,0,0,1;
            0,0,0,0,0,1];
elseif n==3
    a1 = 450e-3;
    a2 = 260e-3;
        % 
    DH = [a1,       0,        0,      0,      0,        1;
          a2,       pi,       0,      0,      0,        1;
          0,        0,        0,      0,      0,        1];
  isRevo=[0,0,0,0,0,1;
          0,0,0,0,0,1;
          0,0,0,0,0,1];  
              
elseif n==4
    a1 = 450e-3;
    a2 = 260e-3;
        % 
    DH = [a1,       0,        0,      0,      0,        1;
          a2,       pi,       0,      0,      0,        1;
          0,        0,        0,      0,      1,        1;
          0,        0,        0,      0,      0,        1];
  isRevo=[0,0,0,0,0,1;
          0,0,0,0,0,1;
          0,0,1,0,0,0;
          0,0,0,0,0,1];  
      
elseif n==5
    d1 = 183e-3;
    d4 = 365e-3;
    d6 = 80e-3;
    a1 = 25e-3;
    a2 = -315e-3;
    a3 = -35e-3;
    DH = [a1,       -pi/2,   d1,     0,      0,      1;
          a2,       0,       0,      0,      0,      1;
          a3,       pi/2,    0,      0,      0,      1;
          0,        -pi/2,   d4,     0,      0,      1;
          0,        0,       d6,     0,      0,      1];
    isRevo=[zeros(5,5),ones(5,1)];
elseif n==6
    d1 = 183e-3;
    d4 = 365e-3;
    d6 = 80e-3;
    a1 = 25e-3;
    a2 = -315e-3;
    a3 = -35e-3;
    DH = [a1,       -pi/2,   d1,     0,      0,      1;
          a2,       0,       0,      0,      0,      1;
          a3,       pi/2,    0,      0,      0,      1;
          0,        -pi/2,   d4,     0,      0,      1;
          0,        pi/2,    0,      0,      0,      1;
          0,        0,       d6,     0,      0,      1];
    isRevo=[zeros(6,5),ones(6,1)];               
elseif n==7
    d_bs = 323.5e-3;
    d_se = 350e-3;
    d_ew = 350e-3;
    d_wt = 275e-3;
        % 
    DH = [0,       -pi/2,    d_bs,     0,      0,        1;
          0,       pi/2,       0,      0,      0,        1;
          0,       -pi/2,    d_se,     0,      0,        1;
          0,       pi/2,      0,       0,      0,        1;
          0,       -pi/2,    d_ew,     0,      0,        1;
          0,       pi/2,        0,     0,      0,        1;
          0,          0,     d_wt,     0,      0,        1];
     isRevo=[zeros(7,5),ones(7,1)];
else
    sprintf('This kind of robot should be config in function config_DH.m')
end



end

