function traceTkJkV = cal_traceV_best_num(T_mat,q,opt)
%cal_traceVfast
% Compute the part V of trace by means of sparse coding
%
% INPUT:
% 	Tmat: syms type.
%   q:joint angle vector
% OUTPUT:
%   traceTkJkV: vector. 
%
% Example:
% 	None
%
% NOTE: 
% 	None
%



% DEFINITION
n=length(q);
traceTkJkV=cell(1,n);
nn=int8(1:16);
T_inx_mat=reshape(nn,[4,4]).';
T_trans_inT_inx=reshape(T_inx_mat,[16,1]);
for i=1:n
%     Jk=fetch_para(i,opt);
    Tk=T_mat{i};
    if isa(Tk.V_fast(T_trans_inT_inx(1),1),'sym')
        V1=kron(Tk.V_fast(T_trans_inT_inx(1),:),Tk.V_fast(1,:))*0;
    else
        V1=kron(zeros(size(Tk.V_fast(T_trans_inT_inx(1),:))),...
                    zeros(size(Tk.V_fast(1,:))));
    end
    len_b=length(Tk.V_fast(1,:));
    
    for j=1:16
        datav=Tk.V_fast(j,:);
       [~,nu]=find(datav);
       inxs{j}=nu;
       datavec{j}=datav(nu);
    end
    ataV1=kron_sparse(inxs,datavec,[T_trans_inT_inx(1),T_trans_inT_inx(2),T_trans_inT_inx(3),T_trans_inT_inx(4)],[1,5,9,13],V1,len_b);
    ataV2=kron_sparse(inxs,datavec,[T_trans_inT_inx(1),T_trans_inT_inx(2),T_trans_inT_inx(3),T_trans_inT_inx(4)],[2,6,10,14],V1,len_b);
    ataV3=kron_sparse(inxs,datavec,[T_trans_inT_inx(1),T_trans_inT_inx(2),T_trans_inT_inx(3),T_trans_inT_inx(4)],[3,7,11,15],V1,len_b);
    ataV4=kron_sparse(inxs,datavec,[T_trans_inT_inx(1),T_trans_inT_inx(2),T_trans_inT_inx(3),T_trans_inT_inx(4)],[4,8,12,16],V1,len_b);
    ataV5=kron_sparse(inxs,datavec,[T_trans_inT_inx(5),T_trans_inT_inx(6),T_trans_inT_inx(7),T_trans_inT_inx(8)],[1,5,9,13],V1,len_b);
    ataV6=kron_sparse(inxs,datavec,[T_trans_inT_inx(5),T_trans_inT_inx(6),T_trans_inT_inx(7),T_trans_inT_inx(8)],[2,6,10,14],V1,len_b);
    ataV7=kron_sparse(inxs,datavec,[T_trans_inT_inx(5),T_trans_inT_inx(6),T_trans_inT_inx(7),T_trans_inT_inx(8)],[3,7,11,15],V1,len_b);
    ataV8=kron_sparse(inxs,datavec,[T_trans_inT_inx(5),T_trans_inT_inx(6),T_trans_inT_inx(7),T_trans_inT_inx(8)],[4,8,12,16],V1,len_b);
    ataV9=kron_sparse(inxs,datavec,[T_trans_inT_inx(9),T_trans_inT_inx(10),T_trans_inT_inx(11),T_trans_inT_inx(12)],[1,5,9,13],V1,len_b);
    ataV10=kron_sparse(inxs,datavec,[T_trans_inT_inx(9),T_trans_inT_inx(10),T_trans_inT_inx(11),T_trans_inT_inx(12)],[2,6,10,14],V1,len_b);
    ataV11=kron_sparse(inxs,datavec,[T_trans_inT_inx(9),T_trans_inT_inx(10),T_trans_inT_inx(11),T_trans_inT_inx(12)],[3,7,11,15],V1,len_b);
    ataV12=kron_sparse(inxs,datavec,[T_trans_inT_inx(9),T_trans_inT_inx(10),T_trans_inT_inx(11),T_trans_inT_inx(12)],[4,8,12,16],V1,len_b);
    ataV13=kron_sparse(inxs,datavec,[T_trans_inT_inx(13),T_trans_inT_inx(14),T_trans_inT_inx(15),T_trans_inT_inx(16)],[1,5,9,13],V1,len_b);
    ataV14=kron_sparse(inxs,datavec,[T_trans_inT_inx(13),T_trans_inT_inx(14),T_trans_inT_inx(15),T_trans_inT_inx(16)],[2,6,10,14],V1,len_b);
    ataV15=kron_sparse(inxs,datavec,[T_trans_inT_inx(13),T_trans_inT_inx(14),T_trans_inT_inx(15),T_trans_inT_inx(16)],[3,7,11,15],V1,len_b);
    ataV16=kron_sparse(inxs,datavec,[T_trans_inT_inx(13),T_trans_inT_inx(14),T_trans_inT_inx(15),T_trans_inT_inx(16)],[4,8,12,16],V1,len_b);

    
  
    temp=zeros(10*n,size(ataV1,2));
%     parameters=[m,rr,ivector];
%     Inxes=[i,i+n:i+n+2,i+4*n:i+4*n+5];
if n==2
    temp(i,:)=ataV16;
    temp(i+n+(i-1)*2,:)=ataV4+ataV13;
    temp(i+4*n+5+(i-1)*5,:)=0.5*(ataV1+ataV6-ataV11);
    traceTkJkV{i}=temp;
else
    temp(i,:)=ataV16;
    temp(i+n+(i-1)*2,:)=ataV4+ataV13;
    temp(i+n+1+(i-1)*2,:)=ataV8+ataV14;
    temp(i+n+2+(i-1)*2,:)=ataV12+ataV15;
    
    temp(i+4*n+(i-1)*5,:)=0.5*(-ataV1+ataV6+ataV11);
    temp(i+4*n+1+(i-1)*5,:)=ataV2+ataV5;
    temp(i+4*n+2+(i-1)*5,:)=ataV3+ataV9;
    temp(i+4*n+3+(i-1)*5,:)=0.5*(ataV1-ataV6+ataV11);
    temp(i+4*n+4+(i-1)*5,:)=ataV7+ataV10;
    temp(i+4*n+5+(i-1)*5,:)=0.5*(ataV1+ataV6-ataV11);
    traceTkJkV{i}=temp;
end
    % 10n*mi
%     traceTkJkV{i}=Jk(1,1)*ataV1+Jk(2,1)*ataV2+Jk(3,1)*ataV3+Jk(4,1)*ataV4+...
%                 Jk(1,2)*ataV5+Jk(2,2)*ataV6+Jk(3,2)*ataV7+Jk(4,2)*ataV8+...
%                 Jk(1,3)*ataV9+Jk(2,3)*ataV10+Jk(3,3)*ataV11+Jk(4,3)*ataV12+...
%                 Jk(1,4)*ataV13+Jk(2,4)*ataV14+Jk(3,4)*ataV15+Jk(4,4)*ataV16;
end
end


function inx = cal_vec2vec(nub,nud,len_d)
%UNTITLED2 
%   nub=NonZeroIndex(vecb),nud=NonZeroIndex(vecd)
%   
    inx=[];
    numb=length(nub);
    for i=1:numb
        inx_i=nud+len_d*(nub(i)-1);
        inx=[inx,inx_i];
    end
end

function V=kron_sparse(inx,datavec,inx_b,inx_d,V,len_d)

    % 
    vec1_inx=cal_vec2vec(inx{inx_b(1)},inx{inx_d(1)},len_d);
    vec_data1=kron(datavec{inx_b(1)},datavec{inx_d(1)});
    
    vec2_inx=cal_vec2vec(inx{inx_b(2)},inx{inx_d(2)},len_d);
    vec_data2=kron(datavec{inx_b(2)},datavec{inx_d(2)});

    vec3_inx=cal_vec2vec(inx{inx_b(3)},inx{inx_d(3)},len_d);
    vec_data3=kron(datavec{inx_b(3)},datavec{inx_d(3)});
    
    vec4_inx=cal_vec2vec(inx{inx_b(4)},inx{inx_d(4)},len_d);
    vec_data4=kron(datavec{inx_b(4)},datavec{inx_d(4)});
    
    
    inter12=intersect(vec1_inx,vec2_inx);
    diffset1=setdiff(vec1_inx,inter12);
    diffset2=setdiff(vec2_inx,inter12);
    V(:,diffset1)=vec_data1(vec1_inx==diffset1);
    V(:,diffset2)=vec_data2(vec2_inx==diffset2);
    if isempty(inter12)~=1
        V(:,inter12)=vec_data1(vec1_inx==inter12)+vec_data2(vec2_inx==inter12);
    end
    % 
    set12=unique([vec1_inx,vec2_inx]);
    inter123=intersect(set12,vec3_inx);
    diffset3=setdiff(vec3_inx,inter123);
    V(:,diffset3)=vec_data3(vec3_inx==diffset3);
    if isempty(inter123)~=1
        V(:,inter123)=V(:,inter123)+vec_data3(vec3_inx==inter123);
    end
    % 
    set123=unique([set12,vec3_inx]);
    inter1234=intersect(set123,vec4_inx);
    diffset4=setdiff(vec4_inx,inter1234);
    V(:,diffset4)=vec_data4(vec4_inx==diffset4);
    if isempty(inter1234)~=1
        V(:,inter1234)=V(:,inter1234)+vec_data4(vec4_inx==inter1234);
    end
end

