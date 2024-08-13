function [NonzerosVec,nu] = NonZerosEle(mats,n)
% reshape the mats to a sparse vector, and return sparse index
% 
mat2vecs=reshape(mats.',[n*n,1]).';
[~,nu]=find(mat2vecs);
NonzerosVec=mat2vecs(nu);
end

