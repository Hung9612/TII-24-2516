function newT= MatchandCalVW_best7_num(amat,vec,Noninx,reinx,rel)
%MatchandCalVW_best7 
% Kronecker product will lead to explosive dimention, so we just use the
% real index of V and W to be the matched index which is encoded inversely
% to seek the value. in a word, we just use the real index to finish the
% operation, then by encoding the matched index into original index to
% attain the value that is exactly needed.
% speed. 
%
% 
%
% Input:
% 	amat: matrix of the first matrixC. 
%   vec: sparse vector of each matrix from matrixC and matrixB
%   Noninx: sparse index for nonzero elements in vec
%   reinx: real inx for the those vector after being processed by function
%   "cal_index_vec2vec.m" and "cal_index_A2vec2I.m"
%   rel:relative index coreponding to reinx, also being processed by the
%   two function.
% Output:
% newT: newT.V and newT.W


newmat=[];
subnewVmat{16,1}=[];
subneVmat_fast{16,1}=[];
inxW=reinx.inx_reducW;
tic
for i=1:4
    temp3=reinx.inx_reducV{i};
    temp4=rel.cur_inx_reducV{i};
%     tic
    for j=1:4
%         tic
        newV_row_fast=[];%只保留不为的avec中与W匹配的元素
        newV_row_inWinx=[];%只保留不为的avec中与W匹配的元素
        row2temp3=cell2mat(temp3(j));%取出V原始的列号（即未简化的）的第j行
        row2temp4=cell2mat(temp4(j));%取当前的列号
%         toc
        
        for k=1:length(row2temp3) %对V的每一行的非零元列号和W的非零元列号进行匹配
%             tic
           inxVij=row2temp3(k);
           
           inxOfV_inW=find(inxW==inxVij);%匹配W和每一行非零元的列号，得到其在W的位置
%            toc
%            tic
           if (inxOfV_inW~=0)% 匹配上了
               
               inx_newVjk=inxOfV_inW;%W为基准的地址，从简化后kron的结果选取Vjk来和W匹配
%                inx_Vjk=row2temp4(k);%取出V有用的relative地址
%                
%                curinx2avec=(inx_Vjk+4-j)/4;%将有用的relative索引映射到avec的relative
%                n_vec=length(reinx.CEGIKM);%计算vec的长度
%                inx_amat=ceil(curinx2avec/n_vec);%计算y=x+n_vec*(nz_amat_col(j)-1)的n_amat_col
%                as=amat(i,inx_amat);
%                if as~=0 % 如果as，即aMat矩阵里的为零，那么直接跳过
%                    curinx2vec=curinx2avec-n_vec*(inx_amat-1);%找出在vec的relative地址
% 
%                    realinx2vec=reinx.CEGIKM(curinx2vec);%将vec的相对地址映射为real地址
%                    inx_vec1=ceil(realinx2vec/16);%由vec的real地址，可以得到上一个vec的real地址
%                    inx_M=realinx2vec-16*(inx_vec1-1);%根据这个地址可以找出M中有用的元素
%                    relative_M=Noninx.M==inx_M;
% 
%                    inx_vec2=ceil(inx_vec1/16);
%                    inx_K=inx_vec1-16*(inx_vec2-1);
%                    relative_K=Noninx.K==inx_K;
% 
%                    inx_vec3=ceil(inx_vec2/16);
%                    inx_I=inx_vec2-16*(inx_vec3-1);
%                    relative_I=Noninx.I==inx_I;
% 
%                    inx_vec4=ceil(inx_vec3/16);
%                    inx_G=inx_vec3-16*(inx_vec4-1);
%                    relative_G=Noninx.G==inx_G;
% 
%                    inx_vec5=ceil(inx_vec4/16);
%                    inx_E=inx_vec4-16*(inx_vec5-1);
%                    relative_E=Noninx.E==inx_E;
% 
%                    inx_vec6=ceil(inx_vec5/16);
%                    inx_C=inx_vec5-16*(inx_vec6-1);
%                    relative_C=Noninx.C==inx_C;
% 
% 
%                    cs=vec.C(relative_C);
%                    es=vec.E(relative_E);
%                    gs=vec.G(relative_G);
%                    is=vec.I(relative_I);
%                    ks=vec.K(relative_K);
%                    ms=vec.M(relative_M);
%                    tempVi_in_avec=as*cs*es*gs*is*ks*ms;
                    
                   tempVi_in_avec=encoder_V_fast(row2temp4,reinx,amat,vec,Noninx,i,j,k);

                   if tempVi_in_avec~=0
                        newV_row_fast=[newV_row_fast,tempVi_in_avec];
                        newmat=[newmat,inx_newVjk];%保存W中有用地址的索引
                        newV_row_inWinx=[newV_row_inWinx,inx_newVjk];
                   end
                   
%                end
               
           end
%            toc
        end
%         tic
        subnewVmat{4*(i-1)+j}=newV_row_inWinx;%每行与W匹配的地址
        subneVmat_fast{4*(i-1)+j}=newV_row_fast;%avec值
%         toc
%         toc
    end
%     toc
end
toc
newmat=unique(newmat); %地址的索引到地址，以及W的值三者是对应的
realinx_W=reinx.inx_reducW(newmat); %获取W有用的地址
newW=zeros(length(newmat),1)*vec.B(1);
% tic
for t=1:length(newmat)
    vecWinx1=ceil(realinx_W(t)/16);
    inx_N=realinx_W(t)-16*(vecWinx1-1);
    relative_N=Noninx.N==inx_N;
    
    vecWinx2=ceil(vecWinx1/16);
    inx_L=vecWinx1-16*(vecWinx2-1);
    relative_L=Noninx.L==inx_L;

    
    vecWinx3=ceil(vecWinx2/16);
    inx_J=vecWinx2-16*(vecWinx3-1);
    relative_J=Noninx.J==inx_J;

    
    vecWinx4=ceil(vecWinx3/16);
    inx_H=vecWinx3-16*(vecWinx4-1);
    relative_H=Noninx.H==inx_H;

    
    vecWinx5=ceil(vecWinx4/16);
    inx_F=vecWinx4-16*(vecWinx5-1);
    relative_F=Noninx.F==inx_F;

    
    vecWinx6=ceil(vecWinx5/16);
    inx_D=vecWinx5-16*(vecWinx6-1);
    relative_D=Noninx.D==inx_D;

    
    vecWinx7=ceil(vecWinx6/16);
    inx_B=vecWinx6-16*(vecWinx7-1);
    relative_B=Noninx.B==inx_B;

    
    tempW1=vec.B(relative_B)*vec.D(relative_D)*vec.F(relative_F)*vec.H(relative_H);
    tempW2=vec.J(relative_J)*vec.L(relative_L)*vec.N(relative_N);
    tempW=tempW1*tempW2;
    newW(t)=tempW;
end
% toc
% 计算V
% newV=zeros(16,length(newmat))*newW(1);
newV=zeros(16,length(newmat));

for i=1:16
    tempVinx=subnewVmat{i};
    tempV_row=subneVmat_fast{i};
%     tempW=T.W(tempVinx);%avec不为零对应的W
    for j=1:length(tempVinx)
        inxinW=newmat==tempVinx(j);%可能会有重复
        newV(i,inxinW)=tempV_row(j);
    end
end

[uniW,inxUni,inverinx]=unique(newW);%地址子集
repetinx=inxUni(inverinx);%返回完整Winx，体现其重复的地址
uniV=newV(:,inxUni)*0;
for i=1:length(inxUni)
    [mu,~]=find(repetinx==inxUni(i));
    reptV=newV(:,mu);
    compressrepeV=sum(reptV,2);
    uniV(:,i)=compressrepeV;
end

newT.W=uniW;
newT.V_fast=uniV;

end

%% local function
function tempVi_in_avec = encoder_V_fast(row2temp4,reinx,amat,vec,Noninx,i,j,k)

inx_Vjk=row2temp4(k);%取出V有用的relative地址

curinx2avec=(inx_Vjk+4-j)/4;%将有用的relative索引映射到avec的relative
n_vec=length(reinx.CEGIKM);%计算vec的长度
inx_amat=ceil(curinx2avec/n_vec);%计算y=x+n_vec*(nz_amat_col(j)-1)的n_amat_col
as=amat(i,inx_amat);

if as~=0 % 如果as，即aMat矩阵里的为零，那么直接跳过
   curinx2vec=curinx2avec-n_vec*(inx_amat-1);%找出在vec的relative地址

   realinx2vec=reinx.CEGIKM(curinx2vec);%将vec的相对地址映射为real地址
   inx_vec1=ceil(realinx2vec/16);%由vec的real地址，可以得到上一个vec的real地址
   inx_M=realinx2vec-16*(inx_vec1-1);%根据这个地址可以找出M中有用的元素
   relative_M=Noninx.M==inx_M;

   inx_vec2=ceil(inx_vec1/16);
   inx_K=inx_vec1-16*(inx_vec2-1);
   relative_K=Noninx.K==inx_K;

   inx_vec3=ceil(inx_vec2/16);
   inx_I=inx_vec2-16*(inx_vec3-1);
   relative_I=Noninx.I==inx_I;

   inx_vec4=ceil(inx_vec3/16);
   inx_G=inx_vec3-16*(inx_vec4-1);
   relative_G=Noninx.G==inx_G;

   inx_vec5=ceil(inx_vec4/16);
   inx_E=inx_vec4-16*(inx_vec5-1);
   relative_E=Noninx.E==inx_E;

   inx_vec6=ceil(inx_vec5/16);
   inx_C=inx_vec5-16*(inx_vec6-1);
   relative_C=Noninx.C==inx_C;


   cs=vec.C(relative_C);
   es=vec.E(relative_E);
   gs=vec.G(relative_G);
   is=vec.I(relative_I);
   ks=vec.K(relative_K);
   ms=vec.M(relative_M);
   tempVi_in_avec=as*cs*es*gs*is*ks*ms;

end

end

