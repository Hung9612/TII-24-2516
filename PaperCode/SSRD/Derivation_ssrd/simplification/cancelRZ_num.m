function [kw,kv] = cancelRZ_num(dij_v,dij_w)

ro=find(sum(abs(dij_v),1));%10n*m
if isempty(ro)==0
    dij_v1=dij_v(:,ro);
    dij_w1=dij_w(:,ro);

    col=find(sum(abs(dij_w1),1));
    if isempty(col)==0
    % col=find(any(dij_w1));
    dij_v2=dij_v1(:,col);
    dij_w2=dij_w1(:,col);
    else
         dij_v2=dij_v1;
         dij_w2=dij_w1;
    end

    [kw1,id,ix]=unique(dij_w2.','rows');
    rx=id(ix);
    if isempty(dij_v2)~=1
        kv1=zeros(size(dij_v2,1),size(id,1));
        for i=1:length(id)
           [r,~]=find(rx==id(i));
           kv1(:,i)=sum(dij_v2(:,r),2);
        end

        col=find(sum(abs(kv1),1));%10n*m

        kv=kv1(:,col);
        kw1=kw1.';
        kw=kw1(:,col);
    else
         kv=dij_v2;
         kw=dij_w2;
    end
else
    kv=dij_v;
    kw=dij_w;
end
end

