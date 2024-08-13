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
%         newV_row=T.W-T.W;%新的V维度大小应当和W匹配，注意循环时释放上一个变量的内存
%         newV_row_fast=T.W*0;
        newV_row_fast=zeros(size(T.W));
        row2temp3=cell2mat(temp3(j));%取出V原始的列号（即未简化的）的第j行
        row2temp4=cell2mat(temp4(j));%取当前的列号
        Inx_deco1=[];
        Inx_Wmat1=[];
        for k=1:length(row2temp3) %对V的每一行的非零元列号和W的非零元列号进行匹配
           inxOfV_inW=find(T.inx_reducW==row2temp3(k));%匹配W和每一行非零元的列号，得到其在W的位置
           if (inxOfV_inW)
               inx_newVjk=inxOfV_inW;%W为基准，从简化后kron的结果选取Vjk来和W匹配
               inx_Vjk=row2temp4(k);
               Inx_avec=(row2temp3(k)+4-j)/4;
               Inx_deco1=[Inx_deco1,Inx_avec];
               Inx_Wmat1=[Inx_Wmat1,row2temp3(k)];
%                newV_row(inx_newVjk) = T.V(4*(i-1)+j,inx_Vjk);%行：4*(i-1)+j，列：inx_Vjk
               newV_row_fast(inx_newVjk)=findinavec(inx_Vjk,T.avec,i,j);
           end

        end
%         subnewVmat=[subnewVmat;newV_row];
        subneVmat_fast=[subneVmat_fast;newV_row_fast];
        Inx_deco=[Inx_deco,Inx_deco1];
        Inx_Wmat=[Inx_Wmat,Inx_Wmat1];
        Inx_Wij{v}=Inx_Wmat1;
        v=v+1;
    end
%     disp(['匹配内循环运行时间:',num2str(toc)])
%     sum(subnewVmat-subneVmat_fast,2)
%     newmat=[newmat;subnewVmat];
    newmat_fast=[newmat_fast;subneVmat_fast];
    Inx_deavec{u}=unique(Inx_deco);
    Inx_W{u}=unique(Inx_Wmat);
    u=u+1;
end
% newT.V=newmat;
% temps=newmat_fast
summVfast=sum(newmat_fast);
[~,j]=find(summVfast);
newT.V_fast=newmat_fast(:,j);
newT.W=T.W(:,j).';%以W为基准
newT.Inxavec=Inx_deavec;
newT.InxWmat=Inx_W;
newT.InxWijmat=Inx_Wij;
% newT.V_fast=newmat_fast;
% newT.W=T.W.';%以W为基准

end

