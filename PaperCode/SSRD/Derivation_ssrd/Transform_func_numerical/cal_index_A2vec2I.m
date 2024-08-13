function [inx,cur_inx] = cal_index_A2vec2I(nu,aIndex,p)
%cal_index_A2vec2I 
% Kronecker product for a matrix, vec vector, and In (not in) identity matrix
% However, the dimention of their result increases soaringly.
% So, we only take the nonzeros elments of them to undertake kronecker.
% Here is calculate their sparse index after kronecker product such that we
% don't care their values, which accelerate the speed.
%
% 
%
% Input:
% 	nu: vec's sparse index which has been computed in function
% 	“cal_index_vec2vec.m”. 
%   aIndex: the sparse index of a mtrix.
%   p: the exponent of gain. here is the exponent of 16,since T matrix of
%   robot is 4 by 4.
%
% Output:
% inx: real index result after the three kronecker product
% cur_inx: relative index change (relative to the number of sparse vector
% or matrix exactly have) .
% 
% Example:
% 	  nu=[1,3,6,7];%nonzero element index of vec vector
%     aIndex.i=[1,1,3,4];%sparse matrix can use fucntion "find"
%    aIndex.j=[2,3,4,5];%
%       p=0;% its size is relavant to the number of C matrix, and B matrix
%   % rerurn
%   inx: cell array.
%   cur_inx:cell array.
%

cur_nu=1:length(nu);
row_index={};

for i=1:4
    nz_col=aIndex.n(aIndex.m==i);
    n=length(nz_col);
    index_col=[];
    curent_indea_col=[];
    
    rates=16*power(16,p);
    for j=1:n
        index_col=[index_col,nu+rates*(nz_col(j)-1)];
        curent_indea_col=[curent_indea_col,cur_nu+length(nu)*(nz_col(j)-1)];
    end
    row_index{i}=index_col;
    cur_row_index{i}=curent_indea_col;
end
cur_aKvec_inx=cur_row_index;

inx={};
cur_inx={};
flag=1;
for i=1:16
    if mod(i,4)==1
        
       temp1=cell2mat(row_index(flag));
       temp0=temp1+3*(temp1-1);
       for  j=1:4
           temp2=temp0+j-1;
           row_i{j}=temp2;
           
           temp3=cell2mat(cur_row_index(flag));
           temp4=temp3+3*(temp3-1)+j-1;
           cur_row_i{j}=temp4;
       end
       
       inx{flag}=row_i;
       cur_inx{flag}=cur_row_i;
       flag=flag+1;
       row_i={};
    end
end


end

