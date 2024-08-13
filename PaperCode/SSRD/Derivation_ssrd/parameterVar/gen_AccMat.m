function out = gen_AccMat(x)
%genV
% Intuitively, we can compute the multiplication of symmetric matrix and a
% column vector easily, but is oppositely hard to work out the
% multiplication of row-wise vectorized D and a vector.
% A direct way is that change thsi vector into the special form:
% x=[x1,x2,x3,x4];
% xgen=[x1,x2,x3,x4,0, 0, 0, 0, 0, 0;
%  0, x1,0, 0, x2,x3,x4,0, 0, 0;
%  0, 0, x1,0, 0, x2,0, x3,x4,0;
%  0, 0, 0, x1,0, 0, x2,0, x3,x4]
% this is what this function does.
% and then Dvec=[n^2 by N ],is calculated from prior steps.
% D*x is changed into :
% DmatX(i,:)=sum(xgen(i,:).'.*Dvec,1);
% Here DmatX can be recognized as the row-wise vectorization of D*x.
%
% Input:
% 	x: must be the column vector .
%
%
% Output:
% xgen=
% [x1,x2,x3,x4,0, 0, 0, 0, 0, 0;
%  0, x1,0, 0, x2,x3,x4,0, 0, 0;
%  0, 0, x1,0, 0, x2,0, x3,x4,0;
%  0, 0, 0, x1,0, 0, x2,0, x3,x4]
%
%

%




n=length(x);
% n0=int(n*(n+1)/2);%int(-1+sqrt(1+8*n));
out=[];
for i =1:n-1
%     out(1,1:n)=x;
%     t=i+1;
%     out(1+i)
aa=zeros(i-1,n-i+1);
bb=x(i:end).';
cc=[zeros(n-i,1),eye(n-i)*x(i)];
    temp=[aa;bb;cc];
    out=[out,temp];
end
temp1=[zeros(n-1,1);x(end)];
out=[out,temp1];

end

