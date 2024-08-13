function [kw,kv] = cancel_nullC_num(cijk_v,cijk_w)

assert(size(cijk_v,1)==size(cijk_w,2),'error')

ro=find(sum(abs(cijk_v),2));
cijk_v1=cijk_v(ro,:);
cijk_w1=cijk_w(:,ro);

col=find(sum(abs(cijk_w1),1));

cijk_v2=cijk_v1(col,:);
cijk_w2=cijk_w1(:,col);


[kw1,id,ix]=unique(cijk_w2.','rows');
rx=id(ix);
% kv1=zeros(length(id),1)*cijk_v2(1,1);
kv1=zeros(length(id),size(cijk_v2,2));

for i=1:length(id)
   [r,~]=find(rx==id(i));
   kv1(i,:)=sum(cijk_v2(r,:),1);

end

col=find(sum(abs(kv1),2));

kv=kv1(col,:).';
kw1=kw1.';
kw=kw1(:,col);


end

