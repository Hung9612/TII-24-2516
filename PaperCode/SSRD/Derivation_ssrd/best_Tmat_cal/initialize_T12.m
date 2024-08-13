function Tmat2 = initialize_T12(Tc,vect,vece)
%UNTITLED10 此处显示有关此函数的摘要
%   此处显示详细说明


[aIndex.m,aIndex.n] = find(Tc{1});%取出a矩阵的非零元的以行序的行号，列号
[V_inx_reduce,~] = cal_index_A2vec2I(vect{2}.inx,aIndex,0);
W_inx_reduce = cal_index_vec2vec(vece{1}.inx,vece{2}.inx);
T.inx_reducV = V_inx_reduce;
T.inx_reducW = W_inx_reduce;
Tmat2 = MatchandCalVW_numT02(T,vect,vece);
end

% local

function Tmat2 = MatchandCalVW_numT02(T,vect,vece)
u=1;
col_W=[];
row_V=[];
for i=1:4
    temp3=T.inx_reducV{i};
    Inx_deco=[];
    Inx_Wmat=[];

    for j=1:4
        row2temp3=cell2mat(temp3(j));%取出V原始的列号（即未简化的）的第j行
        Inx_deco1=[];
        Inx_Wmat1=[];
        newVij=[];
        newWij=[];
        match_inxij = intersect(row2temp3,T.inx_reducW);
        if isempty(match_inxij)~=1
            for k=1:length(match_inxij)
                % 解码V
                matInx_avec = (match_inxij(k)+4-j)/4;%解码到总的avec
                matInx_avec_x = ceil((matInx_avec-0.5)/16); % (matInx_avec_x-1)*16+vec_decode=match_Inxij;
                vec_decode = matInx_avec - (matInx_avec_x-1)*16;
                mval_vec = vect{2}.val(vect{2}.inx == vec_decode)...
                    *vect{1}.val(vect{1}.inx==(matInx_avec_x+(i-1)*4));
                
                % 解码W
                matInx_vec2vec = ceil((match_inxij(k)-0.5)/16);
                vec2_decode = match_inxij(k) - (matInx_vec2vec-1)*16;
                mval_vec2vec = vece{2}.val(vece{2}.inx == vec2_decode)...
                    *vece{1}.val(vece{1}.inx == matInx_vec2vec);
                
                 newVij = [newVij,mval_vec];
                 newWij = [newWij,mval_vec2vec];
                 Inx_deco1=[Inx_deco1,matInx_avec];

            end
        else
            newVij = 0;
            newWij = sym(0);
        end

           u=u+1;
           row_V=blkdiag(row_V,newVij);
           col_W=vertcat(col_W,newWij.');

            Inx_deco=[Inx_deco,Inx_deco1];
            Inx_Wmat=[Inx_Wmat,match_inxij];
    end

     Inx_deavec{i}=unique(Inx_deco);
     Inx_W{i}=unique(Inx_Wmat);
end
[uni_colW,uninx,relinx]=unique(col_W);
reinx=uninx(relinx);
if size(uni_colW,1)~=1&&size(row_V,1)~=1
    row_V_uni = zeros(size(row_V,1),size(uni_colW,1))*row_V(1,1);
end
for tt = 1:length(uninx)
    row_V_uni(:,tt)=sum(row_V(:,reinx==uninx(tt)),2);
end
Tmat2.V_fast = row_V_uni(:,2:end);
Tmat2.W = uni_colW(2:end,:);
Tmat2.Inxavec = Inx_deavec;
Tmat2.InxWmat = Inx_W;
end




