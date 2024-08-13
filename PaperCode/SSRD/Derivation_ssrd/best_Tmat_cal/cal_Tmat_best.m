function Tmat= cal_Tmat_best(theta)

[Tc,Ev] = ConfigTrans(theta);
[vect,vece] = prep(Tc,Ev);

%
[~,Tmat1]=MatMultiply_num({Tc{1}},{Ev{1}},2);
Tmat2 = initialize_T12(Tc,vect,vece);
n=length(theta);

if 2<n
    Inx_deavec = Tmat2.Inxavec;
    Inx_W = Tmat2.InxWmat;
    for t=3:n
        for i=1:4
            temp=[];
            for j=1:length(Inx_deavec{i})
                temp=[temp,(Inx_deavec{i}(j)-1)*16+vect{t}.inx];
            end
            tildeVi{i}=temp;
        end


        for i=1:4
           temp1 = tildeVi{i};
           temp0 = temp1+3*(temp1-1);
           for  j=1:4
               temp2=temp0+j-1;
               InxV_i{j}=temp2;
           end
           InxV{i}=InxV_i;
           InxV_i={};
        end

        for i=1:4
            tildeWij = [];
            for j = 1:length(Inx_W{i})
                tildeWij = [tildeWij,(Inx_W{i}(j)-1)*16+vece{t}.inx];
            end
            tildeWi{i} = tildeWij;
        end

        u=1;
        col_W=[];
        row_V=[];
        for i=1:4
            tildeWii_vector = tildeWi{i};
            Inx_deco = [];
            Inx_Wmat = [];
            for j=1:4
                InxV_ij_vector = InxV{1,i}{1,j};
                match_Inxij = intersect(tildeWii_vector,InxV_ij_vector);
                newVij = [];
                newWij = [];
                Inx_deco1 = [];
                Inx_Wmat1 = [];
                if isempty(match_Inxij)~=1
                   for k = 1:length(match_Inxij) 
                       matInx_avec = (match_Inxij(k)+4-j)/4;
                       [mval_vec,mval_vec2vec] = decode_val(matInx_avec,match_Inxij(k),vect,vece,i,t);
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
                Inx_Wmat=[Inx_Wmat,match_Inxij];
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

        Tmat{t}.V_fast = row_V_uni(:,2:end);
        Tmat{t}.W = uni_colW(2:end,:);
    end
        
end
Tmat{1}.V_fast=Tmat1.V_fast;
Tmat{1}.W=Tmat1.W;
Tmat{2}.V_fast=Tmat2.V_fast;
Tmat{2}.W=Tmat2.W;


