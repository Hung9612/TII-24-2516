function taug = taucg_num_fast(Tmat,tW,q,Gdir)

g=[Gdir(1),Gdir(2),Gdir(3),0].';
Noninx=find(g);
noninx=(Noninx-1)*4+1:4*Noninx-1;

n=length(q);
X=cal_Yd(n);
Gw=[];
Gv=[];
for i=1:n
    GjW={};
    GjV={};
    for p=i:n

        dti=tW{i,p}.difftw;
        temp=zeros(size(X,2),size(dti,2));
        if n==2
            temp(p,:) = -9.81*Tmat{p}.V_fast(Noninx*4,:);
            temp(p+n+(p-1)*2,:) = -9.81*Tmat{p}.V_fast(noninx(1),:);
        else
            temp(p,:) = -9.81*Tmat{p}.V_fast(Noninx*4,:);
            temp(p+n+(p-1)*2:p+n+2+(p-1)*2,:) = -9.81*Tmat{p}.V_fast(noninx,:);            
        end
        GjV{p}=dti(1,:).*temp;
        GjW{p}=dti(2:end,:);
    end

    Gjv=horzcat(GjV{:});
    Gjw=horzcat(GjW{:});

    Gw=[Gw,Gjw];
    sprintf('cal %d Gqi',i)
    Gv=blkdiag(Gv,Gjv.');

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

