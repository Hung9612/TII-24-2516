function [T,newT] = MatMultiply_num(MatrixC,MatrixB,opt)
%MATMULTIPLY 
%  
T=struct;
N=length(MatrixC)+length(MatrixB);
n=length(MatrixB{1});
if opt==1
    coeY=MatrixC{1};
    V=kron(coeY,eye(n));
    baseX=reshape(MatrixB{1}.',[n*n,1]).';
    W=baseX.';
    if N==2 
        T.V=V;
        T.W=W;
    else
        kc=2;
        jb=2;
        for i=2:N-1
            if mod(i,2)==0 
                   kc = kc+1; 
            else 
                
                temp1=reshape(MatrixC{kc-1}.',[n*n,1]).';
                coeY=kron(coeY,temp1);
                V=kron(coeY,eye(n));

                baseX=kron(W.',reshape(MatrixB{jb}.',[n*n,1]).');
                W=baseX.';
                jb=jb+1; 
            end 
            
        end
        T.V=V;
        T.W=W;
    end
else
    [aIndex.m,aIndex.n]=find(MatrixC{1});
    if N==2 
        coeY=MatrixC{1};
        V=kron(coeY,eye(n));
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

        [V_inx_reduce,cur_V_inx]=cal_index_A2vec2I(vec.nu,aIndex,0);
        
        Vecb=reshape(MatrixB{1}.',[n*n,1]).';
        [b.m,b.n]=find(Vecb);
        Vecd=reshape(MatrixB{2}.',[n*n,1]).';
        [d.m,d.n]=find(Vecd);
        W3=kron(Vecb(b.n),Vecd(d.n));
        W_inx_reduce=cal_index_vec2vec(b.n,d.n);
        
        T.W=W3;
        T.inx_reducV=V_inx_reduce;
        T.inx_reducW=W_inx_reduce;
        T.cur_inx_reducV=cur_V_inx;
        T.avec=temp1;
        newT=MatchandCalVW_num1(T);

    elseif N==6
       
        aMat=MatrixC{1};
        
        
        tempC.data=reshape(MatrixC{2}.',[n*n,1]).';
        tempE.data=reshape(MatrixC{3}.',[n*n,1]).';
        [tempC.mu,tempC.nu]=find(tempC.data);
        [tempE.mu,tempE.nu]=find(tempE.data);
        
        % 求简化后的vecc，vecd向量kron
        vecc.data=kron(tempC.data(tempC.nu),tempE.data(tempE.nu));
        vecc.nu=cal_index_vec2vec(tempC.nu,tempE.nu);
       
        
        %计算a,vec,I的kron
        temp1=kron(aMat,vecc.data);
%         V3=kron(temp1,eye(4));% 计算简化后的a和vecc和In的kron积
        [V_inx_reduce,cur_V_inx]=cal_index_A2vec2I(vecc.nu,aIndex,1);
        
        
        % 求vec:b\d\e向量，并简化
        tempB.data=reshape(MatrixB{1}.',[n*n,1]).';
        tempD.data=reshape(MatrixB{2}.',[n*n,1]).'; 
        tempF.data=reshape(MatrixB{3}.',[n*n,1]).';
        [tempB.mu,tempB.nu]=find(tempB.data);
        [tempD.mu,tempD.nu]=find(tempD.data);
        [tempF.mu,tempF.nu]=find(tempF.data);
        
        val1=kron(tempB.data(tempB.nu),tempD.data(tempD.nu));
        val2=kron(val1,tempF.data(tempF.nu));
        W3=val2;% 计算简化后的vecb和vecd,vecf的kron积
        inx1=cal_index_vec2vec(tempB.nu,tempD.nu);
        inx2=cal_index_vec2vec(inx1,tempF.nu);
        W_inx_reduce=inx2;
        
%         T.V=V3;
        T.W=W3;
        T.inx_reducV=V_inx_reduce;
        T.inx_reducW=W_inx_reduce;
        T.cur_inx_reducV=cur_V_inx;
%         T.inx_cur_aKvec=cur_aKvec_inx;
        T.avec=temp1;
        %此时的V3和W3没有对应关系，因此通过W3索引的变换进行匹配
        newT=MatchandCalVW_num1(T);
%         
    elseif N==8
        % 第一个系数矩阵
        aMat=MatrixC{1};
        
        %简化每个vec
        tempC.data=reshape(MatrixC{2}.',[n*n,1]).';
        tempE.data=reshape(MatrixC{3}.',[n*n,1]).';
        tempG.data=reshape(MatrixC{4}.',[n*n,1]).';
        [tempC.mu,tempC.nu]=find(tempC.data);
        [tempE.mu,tempE.nu]=find(tempE.data);
        [tempG.mu,tempG.nu]=find(tempG.data);
        
        % 求简化后的vecc，vecd向量kron
        vecce.data=kron(tempC.data(tempC.nu),tempE.data(tempE.nu));
        vecce.nu=cal_index_vec2vec(tempC.nu,tempE.nu);
        vecceg.data=kron(vecce.data,tempG.data(tempG.nu));
        vecceg.nu=cal_index_vec2vec(vecce.nu,tempG.nu);
        
        %计算a,vec,I的kron
        temp1=kron(aMat,vecceg.data);
%         V4=kron(temp1,eye(4));% 计算简化后的a和vecc和In的kron积
        [V_inx_reduce,cur_V_inx]=cal_index_A2vec2I(vecceg.nu,aIndex,2);
        
        % 求vec:b\d\e向量，并简化
        tempB.data=reshape(MatrixB{1}.',[n*n,1]).';
        tempD.data=reshape(MatrixB{2}.',[n*n,1]).'; 
        tempF.data=reshape(MatrixB{3}.',[n*n,1]).';
        tempH.data=reshape(MatrixB{4}.',[n*n,1]).';
        [tempB.mu,tempB.nu]=find(tempB.data);
        [tempD.mu,tempD.nu]=find(tempD.data);
        [tempF.mu,tempF.nu]=find(tempF.data);
        [tempH.mu,tempH.nu]=find(tempH.data);
        
        val1=kron(tempB.data(tempB.nu),tempD.data(tempD.nu));
        val2=kron(val1,tempF.data(tempF.nu));
        val3=kron(val2,tempH.data(tempH.nu));
        
        W4=val3;% 计算简化后的vecb和vecd,vecf的kron积
        inx1=cal_index_vec2vec(tempB.nu,tempD.nu);
        inx2=cal_index_vec2vec(inx1,tempF.nu);
        inx3=cal_index_vec2vec(inx2,tempH.nu);
        W_inx_reduce=inx3;
        
%         T.V=V4;
        T.W=W4;
        T.inx_reducV=V_inx_reduce;
        T.inx_reducW=W_inx_reduce;
        T.cur_inx_reducV=cur_V_inx;
        T.avec=temp1;
        
        %此时的V3和W3没有对应关系，因此通过W3索引的变换进行匹配
        newT=MatchandCalVW_num(T);
        
    elseif N==10
                % 第一个系数矩阵
        aMat=MatrixC{1};
        
        %简化每个vec
        tempC.data=reshape(MatrixC{2}.',[n*n,1]).';
        tempE.data=reshape(MatrixC{3}.',[n*n,1]).';
        tempG.data=reshape(MatrixC{4}.',[n*n,1]).';
        tempI.data=reshape(MatrixC{5}.',[n*n,1]).';
        [tempC.mu,tempC.nu]=find(tempC.data);
        [tempE.mu,tempE.nu]=find(tempE.data);
        [tempG.mu,tempG.nu]=find(tempG.data);
        [tempI.mu,tempI.nu]=find(tempI.data);
        
        % 求简化后的vecc，vecd向量kron
        
        vecce.data=kron(tempC.data(tempC.nu),tempE.data(tempE.nu));
        vecce.nu=cal_index_vec2vec(tempC.nu,tempE.nu);
        
        vecceg.data=kron(vecce.data,tempG.data(tempG.nu));
        vecceg.nu=cal_index_vec2vec(vecce.nu,tempG.nu);
        
        veccegi.data=kron(vecceg.data,tempI.data(tempI.nu));
        veccegi.nu=cal_index_vec2vec(vecceg.nu,tempI.nu);
        
        %计算a,vec,I的kron
        temp1=kron(aMat,veccegi.data);
%         V5=kron(temp1,eye(4));% 计算简化后的a和vecc和In的kron积
        [V_inx_reduce,cur_V_inx]=cal_index_A2vec2I(veccegi.nu,aIndex,3);
        
        % 求vec:b\d\e向量，并简化
        tempB.data=reshape(MatrixB{1}.',[n*n,1]).';
        tempD.data=reshape(MatrixB{2}.',[n*n,1]).'; 
        tempF.data=reshape(MatrixB{3}.',[n*n,1]).';
        tempH.data=reshape(MatrixB{4}.',[n*n,1]).';
        tempJ.data=reshape(MatrixB{5}.',[n*n,1]).';
        [tempB.mu,tempB.nu]=find(tempB.data);
        [tempD.mu,tempD.nu]=find(tempD.data);
        [tempF.mu,tempF.nu]=find(tempF.data);
        [tempH.mu,tempH.nu]=find(tempH.data);
        [tempJ.mu,tempJ.nu]=find(tempJ.data);
        
        val1=kron(tempB.data(tempB.nu),tempD.data(tempD.nu));
        val2=kron(val1,tempF.data(tempF.nu));
        val3=kron(val2,tempH.data(tempH.nu));
        val4=kron(val3,tempJ.data(tempJ.nu));
        
        W5=val4;% 计算简化后的vecb和vecd,vecf的kron积
        inx1=cal_index_vec2vec(tempB.nu,tempD.nu);
        inx2=cal_index_vec2vec(inx1,tempF.nu);
        inx3=cal_index_vec2vec(inx2,tempH.nu);
        inx4=cal_index_vec2vec(inx3,tempJ.nu);
        W_inx_reduce=inx4;
        
%         T.V=V5;
        T.W=W5;
        T.inx_reducV=V_inx_reduce;
        T.inx_reducW=W_inx_reduce;
        T.cur_inx_reducV=cur_V_inx;
        T.avec=temp1;
        %此时的V3和W3没有对应关系，因此通过W3索引的变换进行匹配
%         newT=MatchandCalVW(T);
        newT=MatchandCalVW_sim(T);
        
    elseif N==12
                        % 第一个系数矩阵
        aMat=MatrixC{1};
        
        %简化每个vec
        tempC.data=reshape(MatrixC{2}.',[n*n,1]).';
        tempE.data=reshape(MatrixC{3}.',[n*n,1]).';
        tempG.data=reshape(MatrixC{4}.',[n*n,1]).';
        tempI.data=reshape(MatrixC{5}.',[n*n,1]).';
        tempK.data=reshape(MatrixC{6}.',[n*n,1]).';
        [tempC.mu,tempC.nu]=find(tempC.data);
        [tempE.mu,tempE.nu]=find(tempE.data);
        [tempG.mu,tempG.nu]=find(tempG.data);
        [tempI.mu,tempI.nu]=find(tempI.data);
        [tempK.mu,tempK.nu]=find(tempK.data);
        
        % 求简化后的vecc，vecd向量kron
        
        vecce.data=kron(tempC.data(tempC.nu),tempE.data(tempE.nu));
        vecce.nu=cal_index_vec2vec(tempC.nu,tempE.nu);
        
        vecceg.data=kron(vecce.data,tempG.data(tempG.nu));
        vecceg.nu=cal_index_vec2vec(vecce.nu,tempG.nu);
        
        veccegi.data=kron(vecceg.data,tempI.data(tempI.nu));
        veccegi.nu=cal_index_vec2vec(vecceg.nu,tempI.nu);
        
        veccegik.data=kron(veccegi.data,tempK.data(tempK.nu));
        veccegik.nu=cal_index_vec2vec(veccegi.nu,tempK.nu);
        
        %计算a,vec,I的kron
        temp1=kron(aMat,veccegik.data);
%         V6=kron(temp1,eye(4));% 计算简化后的a和vecc和In的kron积
        [V_inx_reduce,cur_V_inx]=cal_index_A2vec2I(veccegik.nu,aIndex,4);
        
        % 求vec:b\d\e向量，并简化
        tempB.data=reshape(MatrixB{1}.',[n*n,1]).';
        tempD.data=reshape(MatrixB{2}.',[n*n,1]).'; 
        tempF.data=reshape(MatrixB{3}.',[n*n,1]).';
        tempH.data=reshape(MatrixB{4}.',[n*n,1]).';
        tempJ.data=reshape(MatrixB{5}.',[n*n,1]).';
        tempL.data=reshape(MatrixB{6}.',[n*n,1]).';
        [tempB.mu,tempB.nu]=find(tempB.data);
        [tempD.mu,tempD.nu]=find(tempD.data);
        [tempF.mu,tempF.nu]=find(tempF.data);
        [tempH.mu,tempH.nu]=find(tempH.data);
        [tempJ.mu,tempJ.nu]=find(tempJ.data);
        [tempL.mu,tempL.nu]=find(tempL.data);
        
        val1=kron(tempB.data(tempB.nu),tempD.data(tempD.nu));
        val2=kron(val1,tempF.data(tempF.nu));
        val3=kron(val2,tempH.data(tempH.nu));
        val4=kron(val3,tempJ.data(tempJ.nu));
        val5=kron(val4,tempL.data(tempL.nu));
        
        W6=val5;% 计算简化后的vecb和vecd,vecf的kron积
        inx1=cal_index_vec2vec(tempB.nu,tempD.nu);
        inx2=cal_index_vec2vec(inx1,tempF.nu);
        inx3=cal_index_vec2vec(inx2,tempH.nu);
        inx4=cal_index_vec2vec(inx3,tempJ.nu);
        inx5=cal_index_vec2vec(inx4,tempL.nu);
        W_inx_reduce=inx5;
        
%         T.V=V6;
        T.W=W6;
        T.inx_reducV=V_inx_reduce;
        T.inx_reducW=W_inx_reduce;
        T.cur_inx_reducV=cur_V_inx;
        T.avec=temp1;
        
        %此时的V3和W3没有对应关系，因此通过W3索引的变换进行匹配
%         newT=MatchandCalVW(T);
        newT=MatchandCalVW_sim(T);
    elseif N==14
%         % 第一个系数矩阵
%         aMat=MatrixC{1};
%         %简化每个vec
%         tic
%         tempC.data=reshape(MatrixC{2}.',[n*n,1]).';
%         tempE.data=reshape(MatrixC{3}.',[n*n,1]).';
%         tempG.data=reshape(MatrixC{4}.',[n*n,1]).';
%         tempI.data=reshape(MatrixC{5}.',[n*n,1]).';
%         tempK.data=reshape(MatrixC{6}.',[n*n,1]).';
%         tempM.data=reshape(MatrixC{7}.',[n*n,1]).';
%         [tempC.mu,tempC.nu]=find(tempC.data);
%         [tempE.mu,tempE.nu]=find(tempE.data);
%         [tempG.mu,tempG.nu]=find(tempG.data);
%         [tempI.mu,tempI.nu]=find(tempI.data);
%         [tempK.mu,tempK.nu]=find(tempK.data);
%         [tempM.mu,tempM.nu]=find(tempM.data);
%         toc
%         % 求简化后的vecc，vecd向量kron
%         tic
%         vecce.data=kron(tempC.data(tempC.nu),tempE.data(tempE.nu));
%         vecce.nu=cal_index_vec2vec(tempC.nu,tempE.nu);
%         
%         vecceg.data=kron(vecce.data,tempG.data(tempG.nu));
%         vecceg.nu=cal_index_vec2vec(vecce.nu,tempG.nu);
%         
%         veccegi.data=kron(vecceg.data,tempI.data(tempI.nu));
%         veccegi.nu=cal_index_vec2vec(vecceg.nu,tempI.nu);
%         
%         veccegik.data=kron(veccegi.data,tempK.data(tempK.nu));
%         veccegik.nu=cal_index_vec2vec(veccegi.nu,tempK.nu);
%         
%         veccegikm.data=kron(veccegik.data,tempM.data(tempM.nu));
%         veccegikm.nu=cal_index_vec2vec(veccegik.nu,tempM.nu);
%         toc
%         %计算a,vec,I的kron
% %         tic
% %         temp1=kron(aMat,veccegikm.data);%这里计算时间很长。。。。
% %         for i=1:4
% %             [~,j]=find(aMat(i));
% %             useaMat=aMat(i,j);
% %             for k=length(j)
% %                 
% %             end
% %         end
% %         V7=kron(temp1,eye(4));% 计算简化后的a和vecc和In的kron积
% %         toc
%         [V_inx_reduce,cur_V_inx]=cal_index_A2vec2I(veccegikm.nu,aIndex,4);
%         
%       
%         % 求vec:b\d\e向量，并简化
%         tic
%         tempB.data=reshape(MatrixB{1}.',[n*n,1]).';
%         tempD.data=reshape(MatrixB{2}.',[n*n,1]).'; 
%         tempF.data=reshape(MatrixB{3}.',[n*n,1]).';
%         tempH.data=reshape(MatrixB{4}.',[n*n,1]).';
%         tempJ.data=reshape(MatrixB{5}.',[n*n,1]).';
%         tempL.data=reshape(MatrixB{6}.',[n*n,1]).';
%         tempN.data=reshape(MatrixB{7}.',[n*n,1]).';
%         [tempB.mu,tempB.nu]=find(tempB.data);
%         [tempD.mu,tempD.nu]=find(tempD.data);
%         [tempF.mu,tempF.nu]=find(tempF.data);
%         [tempH.mu,tempH.nu]=find(tempH.data);
%         [tempJ.mu,tempJ.nu]=find(tempJ.data);
%         [tempL.mu,tempL.nu]=find(tempL.data);
%         [tempN.mu,tempN.nu]=find(tempN.data);
%         toc
%         tic
% %         val1=kron(tempB.data(tempB.nu),tempD.data(tempD.nu));
% %         val2=kron(val1,tempF.data(tempF.nu));
% %         val3=kron(val2,tempH.data(tempH.nu));
% %         val4=kron(val3,tempJ.data(tempJ.nu));
% %         val5=kron(val4,tempL.data(tempL.nu));
% %         val6=kron(val5,tempN.data(tempN.nu));
% %         
% %         W7=val6;% 计算简化后的vecb和vecd,vecf的kron积
%         inx1=cal_index_vec2vec(tempB.nu,tempD.nu);
%         inx2=cal_index_vec2vec(inx1,tempF.nu);
%         inx3=cal_index_vec2vec(inx2,tempH.nu);
%         toc
%         tic
%         inx4=cal_index_vec2vec(inx3,tempJ.nu);
%         inx5=cal_index_vec2vec(inx4,tempL.nu);
%         inx6=cal_index_vec2vec(inx5,tempN.nu);
%         W_inx_reduce=inx6;
% %         T.V=V7;
% %         T.W=W7;
%         toc
%         T.inx_reducV=V_inx_reduce;
%         T.inx_reducW=W_inx_reduce;
%         T.cur_inx_reducV=cur_V_inx;
% %         T.avec=temp1;
%         t.n=length(veccegikm.nu);
%         t.vec=veccegikm.data;
%         t.amat=aMat;
        
        %此时的V3和W3没有对应关系，因此通过W3索引的变换进行匹配
%         newT=MatchandCalVW(T);
%          newT=MatchandCalVW_better(T,t);
         
         
         
%%% %%%%%
        %简化每个vec
%         tic
%         tempC.data=reshape(MatrixC{2}.',[n*n,1]).';
%         tempE.data=reshape(MatrixC{3}.',[n*n,1]).';
%         tempG.data=reshape(MatrixC{4}.',[n*n,1]).';
%         tempI.data=reshape(MatrixC{5}.',[n*n,1]).';
%         tempK.data=reshape(MatrixC{6}.',[n*n,1]).';
%         tempM.data=reshape(MatrixC{7}.',[n*n,1]).';
%         [tempC.mu,tempC.nu]=find(tempC.data);
%         [tempE.mu,tempE.nu]=find(tempE.data);
%         [tempG.mu,tempG.nu]=find(tempG.data);
%         [tempI.mu,tempI.nu]=find(tempI.data);
%         [tempK.mu,tempK.nu]=find(tempK.data);
%         [tempM.mu,tempM.nu]=find(tempM.data);
        % 求简化后的vecc，vecd向量kron
%         tic
%         vecce.data=kron(tempC.data(tempC.nu),tempE.data(tempE.nu));
%         vecce.nu=cal_index_vec2vec(tempC.nu,tempE.nu);
        
%         vecceg.data=kron(vecce.data,tempG.data(tempG.nu));
%         vecceg.nu=cal_index_vec2vec(vecce.nu,tempG.nu);
        
%         veccegi.data=kron(vecceg.data,tempI.data(tempI.nu));
%         veccegi.nu=cal_index_vec2vec(vecceg.nu,tempI.nu);
        
%         veccegik.data=kron(veccegi.data,tempK.data(tempK.nu));
%         veccegik.nu=cal_index_vec2vec(veccegi.nu,tempK.nu);
        
%         veccegikm.data=kron(veccegik.data,tempM.data(tempM.nu));
%         veccegikm.nu=cal_index_vec2vec(veccegik.nu,tempM.nu);
%         toc
        %计算a,vec,I的kron
%         tic
%         temp1=kron(aMat,veccegikm.data);%这里计算时间很长。。。。
%         for i=1:4
%             [~,j]=find(aMat(i));
%             useaMat=aMat(i,j);
%             for k=length(j)
%                 
%             end
%         end
%         V7=kron(temp1,eye(4));% 计算简化后的a和vecc和In的kron积
%         toc
        
      
        % 求vec:b\d\e向量，并简化
%         tic
%         tempB.data=reshape(MatrixB{1}.',[n*n,1]).';
%         tempD.data=reshape(MatrixB{2}.',[n*n,1]).'; 
%         tempF.data=reshape(MatrixB{3}.',[n*n,1]).';
%         tempH.data=reshape(MatrixB{4}.',[n*n,1]).';
%         tempJ.data=reshape(MatrixB{5}.',[n*n,1]).';
%         tempL.data=reshape(MatrixB{6}.',[n*n,1]).';
%         tempN.data=reshape(MatrixB{7}.',[n*n,1]).';
%         [tempB.mu,tempB.nu]=find(tempB.data);
%         [tempD.mu,tempD.nu]=find(tempD.data);
%         [tempF.mu,tempF.nu]=find(tempF.data);
%         [tempH.mu,tempH.nu]=find(tempH.data);
%         [tempJ.mu,tempJ.nu]=find(tempJ.data);
%         [tempL.mu,tempL.nu]=find(tempL.data);
%         [tempN.mu,tempN.nu]=find(tempN.data);
%         toc
        
        %计算V的各个vecdata,amat，以及非零元地址real，用于匹配
        aMat=MatrixC{1};
        [vecdata.C,Noninx.C]=NonZerosEle(MatrixC{2},n);
        [vecdata.E,Noninx.E]=NonZerosEle(MatrixC{3},n);
        [vecdata.G,Noninx.G]=NonZerosEle(MatrixC{4},n);
        [vecdata.I,Noninx.I]=NonZerosEle(MatrixC{5},n);
        [vecdata.K,Noninx.K]=NonZerosEle(MatrixC{6},n);
        [vecdata.M,Noninx.M]=NonZerosEle(MatrixC{7},n);
        
        %计算每个vec的real，以及avec的real，relative 地址
        
        realinx.CE=cal_index_vec2vec(Noninx.C,Noninx.E);
        realinx.CEG=cal_index_vec2vec(realinx.CE,Noninx.G);
        realinx.CEGI=cal_index_vec2vec(realinx.CEG,Noninx.I);
        realinx.CEGIK=cal_index_vec2vec(realinx.CEGI,Noninx.K);
        realinx.CEGIKM=cal_index_vec2vec(realinx.CEGIK,Noninx.M);
        [realinx.inx_reducV,relative.cur_inx_reducV]=cal_index_A2vec2I(realinx.CEGIKM,aIndex,5);%这个5很重要
        
   
        %计算W的各个vec，以及非零元地址real
        
        [vecdata.B,Noninx.B]=NonZerosEle(MatrixB{1},n);
        [vecdata.D,Noninx.D]=NonZerosEle(MatrixB{2},n);
        [vecdata.F,Noninx.F]=NonZerosEle(MatrixB{3},n);
        [vecdata.H,Noninx.H]=NonZerosEle(MatrixB{4},n);
        [vecdata.J,Noninx.J]=NonZerosEle(MatrixB{5},n);
        [vecdata.L,Noninx.L]=NonZerosEle(MatrixB{6},n);
        [vecdata.N,Noninx.N]=NonZerosEle(MatrixB{7},n);
        
        %计算W的每个vec的real地址，用于匹配
        realinx.BD=cal_index_vec2vec(Noninx.B,Noninx.D);
        realinx.BDF=cal_index_vec2vec(realinx.BD,Noninx.F);
        realinx.BDFH=cal_index_vec2vec(realinx.BDF,Noninx.H);
        realinx.BDFHJ=cal_index_vec2vec(realinx.BDFH,Noninx.J);
        realinx.BDFHJL=cal_index_vec2vec(realinx.BDFHJ,Noninx.L);
        realinx.inx_reducW=cal_index_vec2vec(realinx.BDFHJL,Noninx.N);
        
        % V: amat,vecdata,Noninx,realinx,relative
        % W: vecdata, Noninx, realinx
        
        newT= MatchandCalVW_best(aMat,vecdata,Noninx,realinx,relative);
        
        
%         tic
%         val1=kron(tempB.data(tempB.nu),tempD.data(tempD.nu));
%         val2=kron(val1,tempF.data(tempF.nu));
%         val3=kron(val2,tempH.data(tempH.nu));
%         val4=kron(val3,tempJ.data(tempJ.nu));
%         val5=kron(val4,tempL.data(tempL.nu));
%         val6=kron(val5,tempN.data(tempN.nu));
%         
%         W7=val6;% 计算简化后的vecb和vecd,vecf的kron积
%         inx1=cal_index_vec2vec(tempB.nu,tempD.nu);
%         inx2=cal_index_vec2vec(inx1,tempF.nu);
%         inx3=cal_index_vec2vec(inx2,tempH.nu);
%         toc
%         tic
%         inx4=cal_index_vec2vec(inx3,tempJ.nu);
%         inx5=cal_index_vec2vec(inx4,tempL.nu);
%         inx6=cal_index_vec2vec(inx5,tempN.nu);
%         W_inx_reduce=inx6;
%         T.V=V7;
%         T.W=W7;
%         toc
%         T.inx_reducV=V_inx_reduce;
%         T.inx_reducW=W_inx_reduce;
%         T.cur_inx_reducV=cur_V_inx;
% %         T.avec=temp1;
%         t.n=length(veccegikm.nu);
%         t.vec=veccegikm.data;
%         t.amat=aMat;
%         
%         %此时的V3和W3没有对应关系，因此通过W3索引的变换进行匹配
% %         newT=MatchandCalVW(T);
%          newT=MatchandCalVW_better(T,t);
    else 
        sprintf("the number of input is error!")
    end

end

end

