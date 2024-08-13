function [Dmat,dtraceW,ddtracew,traceV] = cal_dmat(q,T,opt,pn)
%cal_Dmatnew 
% Dmat is inertia matrix of robot manipulator  
%
%
%
% Input:
% 	q: joint angle vector in syms form,eg. [q1,q2,q3,q4,...]. 
%   T: transform matrix in the pattern of V W
%
%
% Output:
% Dmat: indirect D matrix, that is V for  vector which is made of row-wise
% vectorized D.
% totalV: it is about dynamic parameters, a vector.
% 
% Example:
% 	  None
%   
%
%

% 
n=length(q);
Dmat{n*(1+n)/2,1}=[];


traceV = cal_traceV_best_num(T,q,opt);% 10n*mi

[dtraceW,ddtracew] = cal_traceWnew(T,q,pn);



u=1;


for i=1:n
    for j=i:n
        
        Dmat{u}=cal_dij_num(i,j,n,traceV,dtraceW,pn);% V: m*10n
        if Dmat{u}.W==0
            Dmat{u}.W=zeros(size(Dmat{1}.W,1),1);
        end
        sprintf('cal %d dij',u)
        u=u+1;
    end
end


end

