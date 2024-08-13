function newT= MatchandCalVW_best6_num(amat,vec,Noninx,reinx,rel)


newmat=[];
subnewVmat{16,1}=[];
subneVmat_fast{16,1}=[];
inxW=reinx.inx_reducW;
tic
for i=1:4
    temp3=reinx.inx_reducV{i};
    temp4=rel.cur_inx_reducV{i};
    for j=1:4
        newV_row_fast=[];
        newV_row_inWinx=[];
        row2temp3=cell2mat(temp3(j));
        row2temp4=cell2mat(temp4(j));
        
        for k=1:length(row2temp3) 
            inxVij=row2temp3(k);
           inxOfV_inW=find(inxW==inxVij);
           if (inxOfV_inW~=0)
               
               inx_newVjk=inxOfV_inW;



                    
                   tempVi_in_avec=encoder_V_fast(row2temp4,reinx,amat,vec,Noninx,i,j,k);

                   if tempVi_in_avec~=0
                        newV_row_fast=[newV_row_fast,tempVi_in_avec];
                        newmat=[newmat,inx_newVjk];
                        newV_row_inWinx=[newV_row_inWinx,inx_newVjk];
                   end
                                  
           end
        end

        subnewVmat{4*(i-1)+j}=newV_row_inWinx;
        subneVmat_fast{4*(i-1)+j}=newV_row_fast;

    end
end
toc
newmat=unique(newmat); 
realinx_W=reinx.inx_reducW(newmat); 
newW=zeros(length(newmat),1)*vec.B(1);
for t=1:length(newmat)

    vecWinx2=ceil(realinx_W(t)/16);
    inx_L=realinx_W(t)-16*(vecWinx2-1);
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
    tempW2=vec.J(relative_J)*vec.L(relative_L);%*vec.N(relative_N);
    tempW=tempW1*tempW2;
    newW(t)=tempW;
end

newV=zeros(16,length(newmat));

for i=1:16
    tempVinx=subnewVmat{i};
    tempV_row=subneVmat_fast{i};

    for j=1:length(tempVinx)
        inxinW=newmat==tempVinx(j);
        newV(i,inxinW)=tempV_row(j);
    end
end

[uniW,inxUni,inverinx]=unique(newW);
repetinx=inxUni(inverinx);
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

inx_Vjk=row2temp4(k);

curinx2avec=(inx_Vjk+4-j)/4;
n_vec=length(reinx.CEGIK);

inx_amat=ceil(curinx2avec/n_vec);
as=amat(i,inx_amat);

if as~=0 % 
   curinx2vec=curinx2avec-n_vec*(inx_amat-1);



   realinx2vec=reinx.CEGIK(curinx2vec);
   
   inx_vec2=ceil(realinx2vec/16);
   inx_K=realinx2vec-16*(inx_vec2-1);
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
%    ms=vec.M(relative_M);
%    tempVi_in_avec=as*cs*es*gs*is*ks*ms;
   tempVi_in_avec=as*cs*es*gs*is*ks;


end

end

