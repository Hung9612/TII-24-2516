function checks(twj,tvj,t1,n,q)

base=funcBase(n);
x=cal_Yd(n);
[numX,~] = config_Inertialvalues(n,x);
theta=tvj*x.';
[~,m]=size(twj);
sumw=[];
for i=1:m
    temp=1;
    for j=1:3*n
        e=twj(j,i);
        temp=temp*(base(j))^e;
    end
    sumw=[sumw,temp];
end
err=0;
qnum=rand(5,length(q));
for r=1:2
dijbest=subs(sumw*theta,num2cell([x,q]),num2cell([numX,qnum(r,:)]));
dijt1=subs(t1,num2cell([x,q]),num2cell([numX,qnum(r,:)]));
err=(eval(dijbest)-eval(dijt1))+err;
end
if 1e-9<err
sprintf('------cijk误差为.f%d',err)
end
% simplify(sumw*theta-t1)
end

