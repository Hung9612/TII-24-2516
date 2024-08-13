function [Vn,newV] = filterKronecker(inx,matchedInx,temp1,aKvec_inx,cur_aKvec_inx)


usefulinx={};
for i=1:4
    temp2=inx{i};
    usefulinx1=[];
    for j=1:4
        eachRowInx=cell2mat(temp2(j));
        n=size(eachRowInx,2);
        usefulinx2=zeros(1,n);
        [~,matchj,nonzerosInx]=find(matchedInx(4*(i-1)+j,:));
        for k=1:length(nonzerosInx)
            [~,usefulJ]=find(eachRowInx==nonzerosInx(k));
            usefulinx2(1,usefulJ)=(nonzerosInx(k)-j+1+3)/4;
        end
        usefulinx1=[usefulinx1;usefulinx2];
    end
    usefulinx{i}=usefulinx1;
end
Vn=usefulinx;

cur_inx=[];
new_aKvec_inx={};
for i=1:4
    new_aKvec_inx1=[];
    tempa=Vn{i};
    [~,~,va]=find(tempa);
    vareduce=unique(va);
    num_va=length(vareduce);
    row_inx_aKvec=aKvec_inx{i};
    for j=1:num_va
        [~,matchJ]=find(row_inx_aKvec==vareduce(j));
        cur_inx=[cur_inx,matchJ];
        new_aKvec_inx1=[new_aKvec_inx1,vareduce(j)];
    end
    new_aKvec_inx{i}=new_aKvec_inx1;
end
[uniqueR,uniqueInx,uniqueO]=unique(cur_inx);
new_temp=temp1(:,uniqueR);
newV=kron(new_temp,eye(4));
n=size(new_temp,2);
flag=1;
cur_inx_V={};
for i=1:4
    cur_inxofavec=1:n;
    for j=1:4
        cur_row_inx{j}=cur_inxofavec+3*(cur_inxofavec-1)+j-1;
    end
    cur_inx_V{flag}=cur_row_inx;
    flag=flag+1;
end

end

