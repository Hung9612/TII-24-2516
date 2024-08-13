function tempVi_in_avec = encoder_V_fast_nouse(row2temp4,reinx,amat,vec,Noninx,i,j,k)


inx_Vjk=row2temp4(k);

curinx2avec=(inx_Vjk+4-j)/4;
n_vec=length(reinx.CEGIKM);
inx_amat=ceil(curinx2avec/n_vec);
as=amat(i,inx_amat);

if as~=0 
   curinx2vec=curinx2avec-n_vec*(inx_amat-1);

   realinx2vec=reinx.CEGIKM(curinx2vec);
   inx_vec1=ceil(realinx2vec/16);
   inx_M=realinx2vec-16*(inx_vec1-1);
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

