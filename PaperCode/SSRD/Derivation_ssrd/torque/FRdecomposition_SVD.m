function [B,C] = FRdecomposition_SVD(A,opt)
if opt==4
    [U,S,V]=svd(A);
    r=rank(S);
    sprintf('The number of base parameter is %d',r)
    V=V';
    B=U(:,1:r);
    C=S(1:r,1:r)*V(1:r,:);
elseif opt==7    
    [U,S,V]=svd(A);
    r=rank(S);
        sprintf('The number of base parameter is %d',r)
    V=V';
    B=U(:,1:r);
    C=S(1:r,1:r)*V(1:r,:);
elseif opt==6
    [U,S,V]=svd(A);
    r=rank(S);
        
    V=V';
    B=U(:,1:r)*power(S(1:r,1:r),0.0001);
    C=power(S(1:r,1:r),0.9999)*V(1:r,:);
    er=power(S(1:r,1:r),0.9999)*power(S(1:r,1:r),0.0001)-S(1:r,1:r);
    sprintf('s-er:%.9f',norm(er))

else
   sprintf('error in svd opt') 
end


end

