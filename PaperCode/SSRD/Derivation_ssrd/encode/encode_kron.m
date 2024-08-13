function out = encode_kron(inxw1,inxw2)
%NUMKRON 
[ro,col1]=size(inxw1);
[~,col2]=size(inxw2);

out=zeros(ro,col1*col2);
k=1;
for i=1:col1
    for j=1:col2
        if sum(inxw1(:,i))~=0&&sum(inxw2(:,j))~=0
            out(:,k)=inxw1(:,i)+inxw2(:,j);
        else
            out(:,k)=zeros(ro,1);
        end
        k=k+1;
    end
end


end

