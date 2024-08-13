function taug = taucg_num(Tmat,tW,q,rn)

g=[0,0,-9.81,0].';

Gw=[];
Gv=[];
n=length(q);
X=cal_Yd(n);

for i=1:n
    GjW=[];
    GjV=[];
    tic
    for p=i:n
        [mj,jrc]=cal_Gpara(p,1);
        jrc(4)=mj;
        T=Tmat{p};
        [dti,~,~]=cal_derivite(tW{p}.tw,i,rn,1);
        temp=kron(g.',jrc.')*T.V_fast;
        GjV=[GjV,dti(1,:).*temp];
        GjW=[GjW,dti(2:end,:)];
    end
    Gjv=eval(jacobian(GjV,X));
    toc
    sprintf('cal %d Gqi',i)
    Gw=[Gw,GjW];
    Gv=blkdiag(Gv,Gjv);

end

[uniw,idx,inx]=unique(Gw.','rows');
reinx=idx(inx);
taugv1=zeros(length(idx),size(Gv,2));

for i=1:length(idx)
    [r,~]=find(reinx==idx(i));
    taugv1(i,:)=sum(Gv(r,:),1);
end
taugw1=uniw.';
[ro,~]=find(sum(abs(taugv1),2));

taugv2=taugv1(ro,:);
taugw2=taugw1(:,ro);

[taugw,taugv]=encodesimplify_Vmat(taugv2,taugw2,n);

taug.w=taugw;
taug.v=taugv;

end

function [neww,newv] = encodesimplify_Vmat(kv1,inxw,num)

[univ,idx,inv]=unique(kv1,'rows');
rx=idx(inv);
neww=[];
newv=[];
if num==4
    num=3;
end
for i=1:length(idx)
    [ro,~]=find(rx==idx(i));
    rewinx=inxw(:,ro);
    rrewinx=rewinx;

    for j=1:num

        er1=[2;0]-rrewinx(2*j-1:2*j,:);
        zinx1=find(sum(abs(er1),1)==0);
        
        if zinx1~=0 
            er2=[0;2]-rrewinx(2*j-1:2*j,:);
            zinx2=find(sum(abs(er2),1)==0);
            for t=1:length(zinx1)
                tempp=rrewinx;
                tempp(2*j-1:2*j,:)=tempp(2*j-1:2*j,:)*0;
                err=tempp(:,zinx1(t))-tempp(:,zinx2);
                re=find(sum(abs(err),1)==0);
                if re~=0
                    rrewinx(:,zinx1(t))=tempp(:,zinx1(t));
                    rrewinx(:,zinx2(re))=rrewinx(:,re)*0;

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



