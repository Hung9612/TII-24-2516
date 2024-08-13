function [Cimat,tauc]= taucc_num(dq,Dmat,rn)
%
%   C_i=(q\otimes q)*(Cmat_i)_{vec}
n=length(dq);
Cimat=cell(n*(1+n)/2,n);

tauCW=[];
tauCV=[];
tic
for i=1:n
    u=1;
    ciw=[];
    civ=[];
    for j=1:n
        for k=j:n
        
        Cimat{u,i}=cal_cijkWV_num(Dmat,i,j,k,rn,n);
        
        temp=zeros(n,size(Cimat{u,i}.W,2));
        % 
        if j==k
            temp(j,:)=2;
            tauciw=[Cimat{u,i}.W;temp];
            tauciv=Cimat{u,i}.V;
        else % 
            temp([j,k],:)=1;
            tauciw=[Cimat{u,i}.W;temp];
            tauciv=2*Cimat{u,i}.V;
        end
        ciw=[ciw,tauciw];
%         civ=[civ,tauciv];
        civ=[civ;tauciv];

        u=u+1; 
        end
    end
    [kw,kv]=cancel_nullC_num(civ,ciw);%kv:10n*m
    [ciw1,civ1] = sin2cos2Sim_num(kv,kw,rn);
    tauCW=[tauCW,ciw1];
    tauCV=blkdiag(tauCV,civ1);%m*10nn

end
toc
sprintf('calculate C(q,dq)')
[uniw,idx,inx]=unique(tauCW.','rows');
reinx=idx(inx);

taucv1=zeros(length(idx),size(tauCV,2));

for i=1:length(idx)
    [r,~]=find(reinx==idx(i));
    taucv1(i,:)=sum(tauCV(r,:),1);
end
taucw1=uniw.';
% 
[ro,~]=find(sum(abs(taucv1),2));
taucv2=taucv1(ro,:);
taucw2=taucw1(:,ro);

[taucw,taucv]=encodesimplify_Vmat(taucv2,taucw2,n);

tauc.w=taucw;
tauc.v=taucv;
end

function [neww,newv] = encodesimplify_Vmat(kv1,inxw,n)

[univ,idx,inv]=unique(kv1,'rows');
rx=idx(inv);
neww=[];
newv=[];
if n==4
    n=3;
end
for i=1:length(idx)
    [ro,~]=find(rx==idx(i));
    rewinx=inxw(:,ro);
    rrewinx=rewinx;

    for j=1:n
        % 
        er1=[2;0]-rrewinx(2*j-1:2*j,:);
        zinx1=find(sum(abs(er1),1)==0);%
        % 
        if zinx1~=0 
            er2=[0;2]-rrewinx(2*j-1:2*j,:);
            zinx2=find(sum(abs(er2),1)==0);%
            for t=1:length(zinx1)
                tempp=rrewinx;
                tempp(2*j-1:2*j,:)=tempp(2*j-1:2*j,:)*0;%
                err=tempp(:,zinx1(t))-tempp(:,zinx2);%
                re=find(sum(abs(err),1)==0);
                if re~=0
                    rrewinx(:,zinx1(t))=tempp(:,zinx1(t));
                    rrewinx(:,zinx2(re))=rrewinx(:,re)*0;
%                     inspection(sum(kkw([zinx1(t),zinx2(re)])),rrewinx(:,zinx1(t)));
                end
            end
        end
    end
    nonz=find(sum(rrewinx,1));
    if nonz~=0
        su=0;
    else
        su=1;
        nonz=1;
    end
    neww=[neww,rrewinx(:,nonz)];
    newv=[newv;univ(i,:).*ones(length(nonz),size(univ,2))];

end



end

function inspect_sim(su,nonz,ksw,ro,rrewinx,base,n)
    for m=1:length(nonz)
        temp=1;
        for n=1:2*n
            e=rrewinx(n,nonz(m));
            temp=temp*(base(n))^e;
        end
        su=su+temp;
    end
    er=su-sum(ksw(ro));
    er1=simplify(er);
    sprintf('checkerror：%s',er1)

end


