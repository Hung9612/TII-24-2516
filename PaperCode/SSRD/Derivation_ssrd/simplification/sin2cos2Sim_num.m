function [neww,newv] = sin2cos2Sim_num(kv1,inxw,n)
%ENCODESIMPILFY reduce the terms meeting (sinx)^2+(cosx)^2=1
%   
% [univ,idx,inv]=unique(kv1.','rows');
assert(size(kv1,2)==size(inxw,2),'error')%列为个数

[univ,idx,inv]=unique(kv1.','rows');

rx=idx(inv);
neww=[];
newv=[];
for i=1:length(idx)
    [ro,~]=find(rx==idx(i));
    rewinx=inxw(:,ro);
    rrewinx=rewinx;
%     kkw=ksw(ro);
%     inspection(ksw(ro),rrewinx);
    for j=1:n
        % 先取出相邻两行，即同一个关节角下的sinq, cosq
        er1=[2;0]-rrewinx(2*j-1:2*j,:);
        zinx1=find(sum(abs(er1),1)==0);%取出所有可能潜在可约的sin
        % 如果两个单项式可约，则同一个关节角下，sin，cos不能同在
        if zinx1~=0 
            
            er2=[0;2]-rrewinx(2*j-1:2*j,:);
            zinx2=find(sum(abs(er2),1)==0);%取出所有可能潜在可约的cos
            for t=1:length(zinx1)
                tempp=rrewinx;
                tempp(2*j-1:2*j,:)=tempp(2*j-1:2*j,:)*0;%把需要约去的关节角进行清零
                err=tempp(:,zinx1(t))-tempp(:,zinx2);%一次只拿一个sin去匹配所有可能可约的cos
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
        nonz=1;% 化简了全部都为0，则只需取
    end
    neww=[neww,rrewinx(:,nonz)];
    newv=[newv;univ(i,:).*ones(length(nonz),size(univ,2))];
end
end

function inspect_sim(su,nonz,ksw,ro,rrewinx,n)
base=basefunction(n);

    for m=1:length(nonz)
        temp=1;
        for p=1:2*n
            e=rrewinx(n,nonz(m));
            temp=temp*(base(p))^e;
        end
        su=su+temp;
    end
    er=su-sum(ksw(ro));
    er1=simplify(er);
    sprintf('检验化简后和未编码：%s',er1)

end

