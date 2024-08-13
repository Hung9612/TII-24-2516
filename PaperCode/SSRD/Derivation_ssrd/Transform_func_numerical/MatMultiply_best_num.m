function [T,newT] = MatMultiply_best_num(MatrixC,MatrixB,opt)
%MATMULTIPLY_best
% some coefficient matrices and variable matrices do alternate product
% whose result's vectorized form(row vectorization) can be extracted as 
% product of T.V and T.W.
% 
%
% Input:
% 	MatrixC: coeffcient,in syms type.
%   MatrixB: variable,in syms type.
%   opt:1,means calculation according to the formula
%   opt:2,means only calculte the needed one
% Output:
% newT: a struct type, which  is newT.fast and newT.W to reconstruct
% the serial matrix product. 
% T:only for opt = 1
%
% Example:
% 	  syms z1;syms z2;syms z3; syms z4; syms z5; syms z6;syms z7;
%     syms q1;syms q2;syms q3; syms q4; syms q5; syms q6;syms q7;
%   
%   e1=[zrot(q1),zeros(1,3)';0,0,0,1];
%   e2=[xrot(q2),zeros(1,3)';0,0,0,1];
%   e3=[zrot(q3),zeros(1,3)';0,0,0,1];
%   e4=[xrot(q4),zeros(1,3)';0,0,0,1];
%   e5=[zrot(q5),zeros(1,3)';0,0,0,1];
%   e6=[xrot(q6),zeros(1,3)';0,0,0,1];
%     
%   T6=transl([0,0,z6]);
%   T5=transl([0,0,z5]);
%   T4=transl([0,0,z4]);
%   T3=transl([0,0,z3]);
%   T2=transl([0,0,z2]);
%   T1=transl([0,0,z1]);
%   [~,T06]=MatMultiply_best({T1,T2,T3,T4,T5,T6},{e1,e2,e3,e4,e5,e6},2);
%   % rerurn
%   newT: struct
%
% NOTE: 
% 	1.here are all the symbolic matrix, it is welcome to try the numerical
% 	ones.
%   2.xrot,zrot are self-defined functions, seen in folder:public_func
%   3.the version is fastest one.
%   4. when it is the even times product
% 
% 
%
% 
%
%
%
%   
T=struct;
N=length(MatrixC)+length(MatrixB);
n=length(MatrixB{1});
if opt==1
    coeY=MatrixC{1};
    V=kron(coeY,eye(n));%
    baseX=reshape(MatrixB{1}.',[n*n,1]).';
    W=baseX.';
    if N==2 % 
        T.V=V;
        T.W=W;
    else
        kc=2;
        jb=2;
        for i=2:N-1
            if mod(i,2)==0 % 
                V=kron(eye(n),MatrixC{kc}.')*V;
                   kc = kc+1; % 
            else %
                %
                temp1=reshape(MatrixC{kc-1}.',[n*n,1]).';
                coeY=kron(coeY,temp1);
                V=kron(coeY,eye(n));

                baseX=kron(W.',reshape(MatrixB{jb}.',[n*n,1]).');
                W=baseX.';
                jb=jb+1; % 
            end 
            
        end
        T.V=V;
        T.W=W;
    end
else
    [aIndex.m,aIndex.n]=find(MatrixC{1});%
    if N==2 % 
        coeY=MatrixC{1};
        V=kron(coeY,eye(n));%
        baseX=reshape(MatrixB{1}.',[n*n,1]);
        W=baseX;
        T.V=V;
        T.W=W;
        [i,~]=find(T.W);
        
        newT.V_fast=T.V(:,i);
        newT.W=T.W(i,:);
    elseif N==4
        aMat=MatrixC{1};
        Vecc=reshape(MatrixC{2}.',[n*n,1]).';
        [vec.mu,vec.nu]=find(Vecc);
        temp1=kron(aMat,Vecc(vec.nu));
%         V3=kron(temp1,eye(4));
        [V_inx_reduce,cur_V_inx]=cal_index_A2vec2I(vec.nu,aIndex,0);
        
        Vecb=reshape(MatrixB{1}.',[n*n,1]).';
        [b.m,b.n]=find(Vecb);
        Vecd=reshape(MatrixB{2}.',[n*n,1]).';
        [d.m,d.n]=find(Vecd);
        W3=kron(Vecb(b.n),Vecd(d.n));
        W_inx_reduce=cal_index_vec2vec(b.n,d.n);
        
%         T.V=V3;
        T.W=W3;
        T.inx_reducV=V_inx_reduce;
        T.inx_reducW=W_inx_reduce;
        T.cur_inx_reducV=cur_V_inx;
        T.avec=temp1;
        
        newT=MatchandCalVW_num(T);

    elseif N==6
        %
        aMat=MatrixC{1};
        [vecdata.C,Noninx.C]=NonZerosEle(MatrixC{2},n);
        [vecdata.E,Noninx.E]=NonZerosEle(MatrixC{3},n);
        
        %
        realinx.CE=cal_index_vec2vec(Noninx.C,Noninx.E);
        [realinx.inx_reducV,relative.cur_inx_reducV]=cal_index_A2vec2I(realinx.CE,aIndex,1);%
        
   
        %
        [vecdata.B,Noninx.B]=NonZerosEle(MatrixB{1},n);
        [vecdata.D,Noninx.D]=NonZerosEle(MatrixB{2},n);
        [vecdata.F,Noninx.F]=NonZerosEle(MatrixB{3},n);
        
        
        realinx.BD=cal_index_vec2vec(Noninx.B,Noninx.D);
        realinx.inx_reducW=cal_index_vec2vec(realinx.BD,Noninx.F);
        
        % V: amat,vecdata,Noninx,realinx,relative
        % W: vecdata, Noninx, realinx
        newT= MatchandCalVW_best3_num(aMat,vecdata,Noninx,realinx,relative);
         
    elseif N==8
        %
        aMat=MatrixC{1};
        [vecdata.C,Noninx.C]=NonZerosEle(MatrixC{2},n);
        [vecdata.E,Noninx.E]=NonZerosEle(MatrixC{3},n);
        [vecdata.G,Noninx.G]=NonZerosEle(MatrixC{4},n);
        
        %
        realinx.CE=cal_index_vec2vec(Noninx.C,Noninx.E);
        realinx.CEG=cal_index_vec2vec(realinx.CE,Noninx.G);
        [realinx.inx_reducV,relative.cur_inx_reducV]=cal_index_A2vec2I(realinx.CEG,aIndex,2);
        
   
        %
        [vecdata.B,Noninx.B]=NonZerosEle(MatrixB{1},n);
        [vecdata.D,Noninx.D]=NonZerosEle(MatrixB{2},n);
        [vecdata.F,Noninx.F]=NonZerosEle(MatrixB{3},n);
        [vecdata.H,Noninx.H]=NonZerosEle(MatrixB{4},n);
        
        %
        realinx.BD=cal_index_vec2vec(Noninx.B,Noninx.D);
        realinx.BDF=cal_index_vec2vec(realinx.BD,Noninx.F);
        realinx.inx_reducW=cal_index_vec2vec(realinx.BDF,Noninx.H);
        
        % V: amat,vecdata,Noninx,realinx,relative
        % W: vecdata, Noninx, realinx
        newT= MatchandCalVW_best4_num(aMat,vecdata,Noninx,realinx,relative);
        
    elseif N==10
       %
        aMat=MatrixC{1};
        [vecdata.C,Noninx.C]=NonZerosEle(MatrixC{2},n);
        [vecdata.E,Noninx.E]=NonZerosEle(MatrixC{3},n);
        [vecdata.G,Noninx.G]=NonZerosEle(MatrixC{4},n);
        [vecdata.I,Noninx.I]=NonZerosEle(MatrixC{5},n);
        
        %
        realinx.CE=cal_index_vec2vec(Noninx.C,Noninx.E);
        realinx.CEG=cal_index_vec2vec(realinx.CE,Noninx.G);
        realinx.CEGI=cal_index_vec2vec(realinx.CEG,Noninx.I);
        [realinx.inx_reducV,relative.cur_inx_reducV]=cal_index_A2vec2I(realinx.CEGI,aIndex,3);
        
   
        %
        [vecdata.B,Noninx.B]=NonZerosEle(MatrixB{1},n);
        [vecdata.D,Noninx.D]=NonZerosEle(MatrixB{2},n);
        [vecdata.F,Noninx.F]=NonZerosEle(MatrixB{3},n);
        [vecdata.H,Noninx.H]=NonZerosEle(MatrixB{4},n);
        [vecdata.J,Noninx.J]=NonZerosEle(MatrixB{5},n);
        
        %
        realinx.BD=cal_index_vec2vec(Noninx.B,Noninx.D);
        realinx.BDF=cal_index_vec2vec(realinx.BD,Noninx.F);
        realinx.BDFH=cal_index_vec2vec(realinx.BDF,Noninx.H);
        realinx.inx_reducW=cal_index_vec2vec(realinx.BDFH,Noninx.J);
        
        % V: amat,vecdata,Noninx,realinx,relative
        % W: vecdata, Noninx, realinx
        newT= MatchandCalVW_best5_num(aMat,vecdata,Noninx,realinx,relative);
        
    elseif N==12
        %
        aMat=MatrixC{1};
        [vecdata.C,Noninx.C]=NonZerosEle(MatrixC{2},n);
        [vecdata.E,Noninx.E]=NonZerosEle(MatrixC{3},n);
        [vecdata.G,Noninx.G]=NonZerosEle(MatrixC{4},n);
        [vecdata.I,Noninx.I]=NonZerosEle(MatrixC{5},n);
        [vecdata.K,Noninx.K]=NonZerosEle(MatrixC{6},n);
        
        %
        realinx.CE=cal_index_vec2vec(Noninx.C,Noninx.E);
        realinx.CEG=cal_index_vec2vec(realinx.CE,Noninx.G);
        realinx.CEGI=cal_index_vec2vec(realinx.CEG,Noninx.I);
        realinx.CEGIK=cal_index_vec2vec(realinx.CEGI,Noninx.K);
        [realinx.inx_reducV,relative.cur_inx_reducV]=cal_index_A2vec2I(realinx.CEGIK,aIndex,4);
        
   
        %
        [vecdata.B,Noninx.B]=NonZerosEle(MatrixB{1},n);
        [vecdata.D,Noninx.D]=NonZerosEle(MatrixB{2},n);
        [vecdata.F,Noninx.F]=NonZerosEle(MatrixB{3},n);
        [vecdata.H,Noninx.H]=NonZerosEle(MatrixB{4},n);
        [vecdata.J,Noninx.J]=NonZerosEle(MatrixB{5},n);
        [vecdata.L,Noninx.L]=NonZerosEle(MatrixB{6},n);
        
        %
        realinx.BD=cal_index_vec2vec(Noninx.B,Noninx.D);
        realinx.BDF=cal_index_vec2vec(realinx.BD,Noninx.F);
        realinx.BDFH=cal_index_vec2vec(realinx.BDF,Noninx.H);
        realinx.BDFHJ=cal_index_vec2vec(realinx.BDFH,Noninx.J);
        realinx.inx_reducW=cal_index_vec2vec(realinx.BDFHJ,Noninx.L);
        
        % V: amat,vecdata,Noninx,realinx,relative
        % W: vecdata, Noninx, realinx
        newT= MatchandCalVW_best6_num(aMat,vecdata,Noninx,realinx,relative);
        
    elseif N==14
        %
        aMat=MatrixC{1};
        [vecdata.C,Noninx.C]=NonZerosEle(MatrixC{2},n);
        [vecdata.E,Noninx.E]=NonZerosEle(MatrixC{3},n);
        [vecdata.G,Noninx.G]=NonZerosEle(MatrixC{4},n);
        [vecdata.I,Noninx.I]=NonZerosEle(MatrixC{5},n);
        [vecdata.K,Noninx.K]=NonZerosEle(MatrixC{6},n);
        [vecdata.M,Noninx.M]=NonZerosEle(MatrixC{7},n);
        
        %
        realinx.CE=cal_index_vec2vec(Noninx.C,Noninx.E);
        realinx.CEG=cal_index_vec2vec(realinx.CE,Noninx.G);
        realinx.CEGI=cal_index_vec2vec(realinx.CEG,Noninx.I);
        realinx.CEGIK=cal_index_vec2vec(realinx.CEGI,Noninx.K);
        realinx.CEGIKM=cal_index_vec2vec(realinx.CEGIK,Noninx.M);
        [realinx.inx_reducV,relative.cur_inx_reducV]=cal_index_A2vec2I(realinx.CEGIKM,aIndex,5);%
        
   
        %
        [vecdata.B,Noninx.B]=NonZerosEle(MatrixB{1},n);
        [vecdata.D,Noninx.D]=NonZerosEle(MatrixB{2},n);
        [vecdata.F,Noninx.F]=NonZerosEle(MatrixB{3},n);
        [vecdata.H,Noninx.H]=NonZerosEle(MatrixB{4},n);
        [vecdata.J,Noninx.J]=NonZerosEle(MatrixB{5},n);
        [vecdata.L,Noninx.L]=NonZerosEle(MatrixB{6},n);
        [vecdata.N,Noninx.N]=NonZerosEle(MatrixB{7},n);
        
        %
        realinx.BD=cal_index_vec2vec(Noninx.B,Noninx.D);
        realinx.BDF=cal_index_vec2vec(realinx.BD,Noninx.F);
        realinx.BDFH=cal_index_vec2vec(realinx.BDF,Noninx.H);
        realinx.BDFHJ=cal_index_vec2vec(realinx.BDFH,Noninx.J);
        realinx.BDFHJL=cal_index_vec2vec(realinx.BDFHJ,Noninx.L);
        realinx.inx_reducW=cal_index_vec2vec(realinx.BDFHJL,Noninx.N);
        
        % V: amat,vecdata,Noninx,realinx,relative
        % W: vecdata, Noninx, realinx
        newT= MatchandCalVW_best7_num(aMat,vecdata,Noninx,realinx,relative);
    else 
        sprintf("the number of input is error!")
    end

end

end



