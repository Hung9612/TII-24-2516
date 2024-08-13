function inx = cal_index_vec2vec(nub,nud)

inx=[];
for i=1:length(nub)
    inx_i=nud+16*(nub(i)-1);
    inx=[inx,inx_i];
end
end

