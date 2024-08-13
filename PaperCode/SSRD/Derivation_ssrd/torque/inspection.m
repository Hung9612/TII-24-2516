function  inspection(diff_tempW,twj,n)

base=funcBase(n);
for i=1:length(diff_tempW)
    temp=1;
    for j=1:3*n
        e=twj(j,i);
        temp=temp*(base(j))^e;
    end
    if temp~=1
       er=simplify(diff_tempW(i)-temp);
       sprintf('%s',er)
    end
end
end

