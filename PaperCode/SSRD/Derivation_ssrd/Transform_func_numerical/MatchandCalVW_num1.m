function newT = MatchandCalVW_num1(T)

newmat_fast=[];
% newmat=[];
u=1;
v=1;
for i=1:4
%     subnewVmat=[];
    subneVmat_fast=[];
    temp3=T.inx_reducV{i};
    temp4=T.cur_inx_reducV{i};
    Inx_deco=[];
    Inx_Wmat=[];

    for j=1:4
%         newV_row=T.W-T.W;
%         newV_row_fast=T.W*0;
        newV_row_fast=zeros(size(T.W));
        row2temp3=cell2mat(temp3(j));
        row2temp4=cell2mat(temp4(j));
        Inx_deco1=[];
        Inx_Wmat1=[];
        for k=1:length(row2temp3) 
           inxOfV_inW=find(T.inx_reducW==row2temp3(k));
           if (inxOfV_inW)
               inx_newVjk=inxOfV_inW;
               inx_Vjk=row2temp4(k);
               Inx_avec=(row2temp3(k)+4-j)/4;
               Inx_deco1=[Inx_deco1,Inx_avec];
               Inx_Wmat1=[Inx_Wmat1,row2temp3(k)];
               newV_row_fast(inx_newVjk)=findinavec(inx_Vjk,T.avec,i,j);
           end

        end
        subneVmat_fast=[subneVmat_fast;newV_row_fast];
        Inx_deco=[Inx_deco,Inx_deco1];
        Inx_Wmat=[Inx_Wmat,Inx_Wmat1];
        Inx_Wij{v}=Inx_Wmat1;
        v=v+1;
    end
    newmat_fast=[newmat_fast;subneVmat_fast];
    Inx_deavec{u}=unique(Inx_deco);
    Inx_W{u}=unique(Inx_Wmat);
    u=u+1;
end

summVfast=sum(newmat_fast);
[~,j]=find(summVfast);
newT.V_fast=newmat_fast(:,j);
newT.W=T.W(:,j).';%
newT.Inxavec=Inx_deavec;
newT.InxWmat=Inx_W;
newT.InxWijmat=Inx_Wij;
% newT.V_fast=newmat_fast;
% newT.W=T.W.';

end

