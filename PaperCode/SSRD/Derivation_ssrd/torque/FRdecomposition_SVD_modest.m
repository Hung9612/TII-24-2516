function [B,C] = FRdecomposition_SVD_modest(Ain,opt)
A=Ain.';
[~,~,V]=svd(A);

if ismember(opt,[2,4,6,7])
    if opt==7&&10000<size(A,1)||10000<size(A,2)
    [R,jb]=rref(V.',1e-6);
    else
    [R,jb]=rref(V.');
    end
    C=A(:,jb).';
    B=R(1:rank(A),:).';
    er=B*C-Ain;
    sprintf('s-er:%.9f',norm(er))
else
   sprintf('error in svd opt') 
end
% if ismember(opt,[2,4,6,7])
%     if opt==7&&10000<size(A,1)||10000<size(A,2)
%     [R,jb]=rref(A,1e-6);
%     else
%     [R,jb]=rref(A);
%     end
%     C=A(:,jb).';
%     B=R(1:rank(A),:).';
%     er=B*C-Ain;
%     sprintf('s-er:%.9f',norm(er))
% else
%    sprintf('error in svd opt') 
% end


end

