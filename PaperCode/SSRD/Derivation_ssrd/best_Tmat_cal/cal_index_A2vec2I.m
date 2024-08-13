function [inx,cur_inx] = cal_index_A2vec2I(nu,aIndex,p)



cur_nu=1:length(nu);
row_index={};

%
for i=1:4
    nz_col=aIndex.n(aIndex.m==i);%
    n=length(nz_col);
    index_col=[];
    curent_indea_col=[];
    
    rates=16*power(16,p);
    for j=1:n
        index_col=[index_col,nu+rates*(nz_col(j)-1)];
        curent_indea_col=[curent_indea_col,cur_nu+length(nu)*(nz_col(j)-1)];
    end
    row_index{i}=index_col;%
    cur_row_index{i}=curent_indea_col;%
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
           row_i{j}=temp2;% 
           
           temp3=cell2mat(cur_row_index(flag));
           temp4=temp3+3*(temp3-1)+j-1;
           cur_row_i{j}=temp4;% 
       end
       
       inx{flag}=row_i;
       cur_inx{flag}=cur_row_i;
       flag=flag+1;
       row_i={};
    end
end


end

