function [B,C] = FRdecomposition_SVD_thin(Ain,opt)
A=Ain.';
if ismember(opt,[2,3,4,5,6,7])
    if opt==7&&10000<size(A,1)||10000<size(A,2)
    [R,jb]=rref(A,1e-6);
    else
    [R,jb]=rref(A);
    end
    C=A(:,jb).';
    B=R(1:rank(A),:).';
    er=B*C-Ain;
    sprintf('s-er:%.9f',norm(er))
else
   sprintf('error in svd opt') 
end


end

