function tw = cal_derivite_best(tw,qi,pn)
%CAL_DERIVATE 在编码中求导
%   -tw:变换矩阵的W的编码矩阵[sin(q);cos(q);q]
%   -pn:移动关节的编号
%   -qi:要求导的关节编号
%   -deg:the function degree, namely the power of the function sin and cos
assert(isempty(tw)~=1,'The tw is empty!')
n=(size(tw,1)-1)/3;
if pn==0 % 只有revolute
    sinx=find(tw(2*qi,:));
    cosx=find(tw(2*qi+1,:));

    scinx=union(sinx,cosx);

    comple=1:size(tw,2);
    % remove the term whose partial diff is zero 
    tw(:,setdiff(comple,scinx))=0;

    % composite function differential
    boinx=intersect(sinx,cosx);
    cx=setdiff(cosx,boinx);
    sx=setdiff(sinx,boinx);
    % both have sine and cosine
    if ~isempty(boinx)==1&&~isempty(cosx)==1&&isempty(sinx)~=1
        tw(1,boinx)=tw(1,boinx)*-1;
    end
    % diff sin
    if ~isempty(sx)==1
        tw(2*qi+1,sx)=1;
        tw(2*qi,sx)=0;
    end
    % diff cos
    if ~isempty(cx)==1
        tw(2*qi,cx)=1;
        tw(2*qi+1,cx)=0;
        tw(1,cx)=tw(1,cx)*-1;
    end 
elseif pn~=0&&pn<=n % revolute and prismatrical
    if qi~=pn
        sinx=find(tw(2*qi,:));
        cosx=find(tw(2*qi+1,:));

        scinx=union(sinx,cosx);
        cop=1:size(tw,2);
        % 
        tw(:,setdiff(cop,scinx))=0;

        boinx=intersect(sinx,cosx);
        cx=setdiff(cosx,boinx);
        sx=setdiff(sinx,boinx);
        % 
        if ~isempty(boinx)==1&&~isempty(cosx)==1&&isempty(sinx)~=1
            tw(1,boinx)=tw(1,boinx)*-1;
        end
        if isempty(sx)~=1
            tw(2*qi+1,sx)=1;
            tw(2*qi,sx)=0;
        end
        if isempty(cx)~=1
            tw(2*qi,cx)=1;
            tw(2*qi+1,cx)=0;
            tw(1,cx)=tw(1,cx)*-1;
        end
    else
        qx=find(tw(2*n+1+pn,:));
        if isempty(qx)~=1
            comple=1:size(tw,2);
            %
            tw(:,setdiff(comple,qx))=0;
            tw(2*n+1+pn,qx)=tw(2*n+1+pn,qx)-1;
        else
            tw(1,:)=0;
        end
    end
else
    sprintf('Wrong! May be that pn is greater than n.')
end
end





