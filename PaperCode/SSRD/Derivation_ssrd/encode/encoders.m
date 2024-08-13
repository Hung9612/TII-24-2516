function inxW = encoders(Wfunc,base,rj,pn)
%ENCODENUM encode the function base in Transformation matrix.
%   [sinq;cosq;q]->3*n
m=length(Wfunc);
mfunc=length(base);
n=uint8(mfunc/3);
inxW=zeros(mfunc,m);

if pn==0
    for i=1:m
        temp=Wfunc(i);
        for j=1:mfunc-n
            flag=subs(temp,{base(j)},{0});
            if flag==0&&temp~=0
                inxW(j,i)=1;
            end
        end
    end
elseif pn~=0
    % baseinx=1:3*n;% 把pn旋转相关的函数和平移的函数取出来
    omitinx=[2*pn-1:2*pn,setdiff(2*n+1:3*n,2*n+pn)];

    % base=base(setdiff(baseinx,omitinx));% 默认有平移就没有旋转
    base(omitinx)=0;
    mfunc=length(base);
    for i=1:m
        temp=Wfunc(i);
        for j=1:mfunc
            flag=subs(temp,{base(j)},{0});
            if flag==0&&temp~=0
                inxW(j,i)=1;
            end
        end
    end
end
inxW=[rj;inxW];


end

