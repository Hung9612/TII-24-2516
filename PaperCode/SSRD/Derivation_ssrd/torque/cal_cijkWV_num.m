function Cijk =cal_cijkWV_num(Dmat,i,j,k,rn,n)

tij_k=calderivite([ones(1,size(Dmat{cal_inx(i,j,n)}.W,2));Dmat{cal_inx(i,j,n)}.W],k,rn);
tik_j=calderivite([ones(1,size(Dmat{cal_inx(i,k,n)}.W,2));Dmat{cal_inx(i,k,n)}.W],j,rn);
tjk_i=calderivite([ones(1,size(Dmat{cal_inx(j,k,n)}.W,2));Dmat{cal_inx(j,k,n)}.W],i,rn);

cijk_w=[tjk_i(2:end,:),tik_j(2:end,:),tij_k(2:end,:)];

cjk=-1/2*double(tjk_i(1,:).').*(Dmat{cal_inx(j,k,n)}.V);%行为个数
cik=1/2*double(tik_j(1,:).').*(Dmat{cal_inx(i,k,n)}.V);
cij=1/2*double(tij_k(1,:).').*(Dmat{cal_inx(i,j,n)}.V);
cijk_v=[cjk;cik;cij];
Cijk.W=cijk_w;
Cijk.V=cijk_v;
end


function inx = cal_inx(s,t,n)

inx=(s-1)*(n-s/2)+t;
end


