function newT = MatchandCalVW_num(T)

newmat_fast=[];
% newmat=[];

for i=1:4
%     subnewVmat=[];
    subneVmat_fast=[];
    temp3=T.inx_reducV{i};
    temp4=T.cur_inx_reducV{i};
    for j=1:4
%         newV_row=T.W-T.W;
%         newV_row_fast=T.W*0;
        newV_row_fast=zeros(size(T.W));
        row2temp3=cell2mat(temp3(j));
        row2temp4=cell2mat(temp4(j));
        for k=1:length(row2temp3) 
           inxOfV_inW=find(T.inx_reducW==row2temp3(k));
           if (inxOfV_inW)
               inx_newVjk=inxOfV_inW;
               inx_Vjk=row2temp4(k);
%                newV_row(inx_newVjk) = T.V(4*(i-1)+j,inx_Vjk);
               newV_row_fast(inx_newVjk)=findinavec(inx_Vjk,T.avec,i,j);
           end

        end
       
%         subnewVmat=[subnewVmat;newV_row];
        subneVmat_fast=[subneVmat_fast;newV_row_fast];
    end

    newmat_fast=[newmat_fast;subneVmat_fast];

end
% newT.V=newmat;
% temps=newmat_fast
summVfast=sum(newmat_fast);
[~,j]=find(summVfast);
newT.V_fast=newmat_fast(:,j);
newT.W=T.W(:,j).';%


end

