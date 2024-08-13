function [mval_vec,mval_vec2vec] = decode_val(matInx_avec,match_Inxij,vect,vece,i,t)


mval_vec=1;
mval_vec2vec=1;
matInx_vec2vec = match_Inxij;
for kk=1:t-1
    

    matInx_avec_x = ceil((matInx_avec-0.5)/16); 
    vec_decode = matInx_avec - (matInx_avec_x-1)*16;
    %
    mval_vec = mval_vec * vect{t-kk+1}.val(vect{t-kk+1}.inx == vec_decode);
     %
    matInx_avec = matInx_avec_x;
    

    matInx_vec2vec = ceil((match_Inxij-0.5)/16);
    vec2_decode = match_Inxij - (matInx_vec2vec-1)*16;
    mval_vec2vec = mval_vec2vec * vece{t-kk+1}.val(vece{t-kk+1}.inx == vec2_decode);
    match_Inxij = matInx_vec2vec;
end
mval_vec = mval_vec * vect{1}.val(vect{1}.inx==(matInx_avec_x+(i-1)*4));

mval_vec2vec = mval_vec2vec * vece{1}.val(vece{1}.inx == matInx_vec2vec);
end

