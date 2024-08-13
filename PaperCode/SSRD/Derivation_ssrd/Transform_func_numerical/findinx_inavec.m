function out = findinx_inavec(usecurinx,t,i,j)

curinx2avec=(usecurinx+4-j)/4;
inx_amat=ceil(curinx2avec/t.n);
curinx2vec=curinx2avec-t.n*(inx_amat-1);

out=t.vec(curinx2vec)*t.amat(i,inx_amat);

end

