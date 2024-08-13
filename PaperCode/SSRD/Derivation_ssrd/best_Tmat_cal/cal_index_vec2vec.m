function inx = cal_index_vec2vec(nub,nud)
%cal_index_vec2vec 
% Kronecker product for vector and vector.
% However, the dimention of their result increases soaringly.
% So, we only take the nonzeros elments of them to undertake kronecker.
% Here is calculate their sparse index after kronecker product such that we
% don't care their values, which accelerate the speed.
%
% 中文：计算 kron(vec(b),vec(d))的索引变换
%
% Input:
% 	nub: b vector sparse index. 
%   nud: d vector sparse index
%
% Output:
% inx: real index result after the three kronecker product

% 
% Example:
% 	  nub=[1,3,6,7];%nonzero element index of vecb vector
%    nud=[4,10,67,47];%nonzero element index of vec vector
%   % rerurn
%   inx: vector index.
%   
%
%

%   nub=NonZeroIndex(vecb),nud=NonZeroIndex(vecd)
%   计算 kron(vec(b),vec(d))的索引变换
inx=[];
for i=1:length(nub)
    inx_i=nud+16*(nub(i)-1);
    inx=[inx,inx_i];
end
end

