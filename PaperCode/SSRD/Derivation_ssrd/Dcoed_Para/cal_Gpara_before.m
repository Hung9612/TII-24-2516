function [mi,ri] = cal_Gpara_before(n,opt)
%fetchGpara 
% For this D matrix, it will compute the result with respect to the pseudo 
% inertia matix   
%
%
% 
%
% Input:
%  n:the number of joints
%
% Output:
% mi: is a  scale
% 
% ri: it is a vector.
% 
% Example:
% 	  None
%   
%
%
%
if opt==2
    syms m1;syms m2;syms m3;
    syms m1rx1;syms m1ry1;syms m1rz1;
    syms m2rx2;syms m2ry2;syms m2rz2;
    syms m3rx3;syms m3ry3;syms m3rz3;

    m=[m1,m2,m3];
    r=[0,  0,  m1rz1;
       0,  0,  m2rz2;
       0,  0,  m3rz3];
    mi=m(n);
    ri=r(n,:).';
    ri=[ri;mi];
else
    syms m1;syms m2;syms m3;syms m4;syms m5;syms m6;syms m7;
    syms m1rx1;syms m1ry1;syms m1rz1;
    syms m2rx2;syms m2ry2;syms m2rz2;
    syms m3rx3;syms m3ry3;syms m3rz3;
    syms m4rx4;syms m4ry4;syms m4rz4;
    syms m5rx5;syms m5ry5;syms m5rz5;
    syms m6rx6;syms m6ry6;syms m6rz6;
    syms m7rx7;syms m7ry7;syms m7rz7;

    m=[m1,m2,m3,m4,m5,m6,m7];
    r=[m1rx1, m1ry1,  m1rz1;
      m2rx2,  m2ry2,  m2rz2;
      m3rx3,  m3ry3,  m3rz3;
      m4rx4,  m4ry4,  m4rz4;
      m5rx5,  m5ry5,  m5rz5;
      m6rx6,  m6ry6,  m6rz6;
      m7rx7,  m7ry7,  m7rz7;];

    %
    mi=m(n);
    ri=r(n,:).';
    ri=[ri;mi];
end
end

