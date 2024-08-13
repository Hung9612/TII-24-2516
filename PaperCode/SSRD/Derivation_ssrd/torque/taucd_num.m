function taud = taucd_num(dmat,ddq)

num=length(ddq);
Mddq=gen_AccMat(ddq.');
taudw=[];
taudv=[];
for i=1:num
    ainx=find(Mddq(i,:));

    tdiw=cell(length(ainx),1);
    tdiv=cell(length(ainx),1);
    for j=1:length(ainx)
        dij=dmat{ainx(j)};
        temp=zeros(num,size(dij.V,1));
        temp(j,:)=1;
        if size(dij.W,2)==size(dij.V,1)% dij.V:m*10n

            tdiw{j}=[dij.W;temp];
            tdiv{j}=dij.V;
        else
            sprintf('error: dimension')
            return
        end
    end
    tdiw1=horzcat(tdiw{:});
    tdiv1=vertcat(tdiv{:});

    taudw=[taudw,tdiw1];
    taudv=blkdiag(taudv,tdiv1);% m*10nn
end

[uniw,idx,inx]=unique(taudw.','rows');
reinx=idx(inx);
% taudv1=zeros(length(idx),num);
taudv1=zeros(length(idx),size(taudv,2));

for i=1:length(idx)
    [r,~]=find(reinx==idx(i));
    taudv1(i,:)=sum(taudv(r,:),1);
end
taudw1=uniw.';

[ro,~]=find(sum(abs(taudv1),2));

% taudv2=taudv1(ro,:);
taudv2=taudv1(ro,:);
taudw2=taudw1(:,ro);

[taudw,taudv]=encodesimplify_Vmat(taudv2,taudw2,num);
taud.w=taudw;
taud.v=taudv;

end

function [neww,newv] = encodesimplify_Vmat(kv1,inxw,num)

if num==4
    num=3;
end
[univ,idx,inv]=unique(kv1,'rows');
rx=idx(inv);
% base=basefunction(7);
neww=[];
newv=[];
for i=1:length(idx)
    [ro,~]=find(rx==idx(i));
    rewinx=inxw(:,ro);
    rrewinx=rewinx;
%     kkw=ksw(ro);
%     inspection(ksw(ro),rrewinx);
    for j=1:num
        % 先取出相邻两行，即同一个关节角下的sinq, cosq
        er1=[2;0]-rrewinx(2*j-1:2*j,:);
        zinx1=find(sum(abs(er1),1)==0);
        % 
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
        su=0;
        nonz=1;
    end
    neww=[neww,rrewinx(:,nonz)];
    newv=[newv;univ(i,:).*ones(length(nonz),size(univ,2))];
end

end

