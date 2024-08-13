function [tw,attachterm,boinx] = cal_derivite(tw,qi,pn,deg)
%CAL_DERIVATE 在编码中求导
%   -tw:变换矩阵的W的编码矩阵[sin(q);cos(q);q]
%   -pn:移动关节的编号
%   -qi:要求导的关节编号
%   -deg:the function degree, namely the power of the function sin and cos
assert(isempty(tw)~=1,'The tw is empty!')
n=(size(tw,1)-1)/3;
attachterm=[];
if deg==1 % cal D(q)
    if pn==0 % 只有revolute
        sinx=find(tw(2*qi,:));
        cosx=find(tw(2*qi+1,:));

        scinx=union(sinx,cosx);

        comple=1:size(tw,2);
        tw(:,setdiff(comple,scinx))=0;

        boinx=intersect(sinx,cosx);
        cx=setdiff(cosx,boinx);
        sx=setdiff(sinx,boinx);
        if ~isempty(boinx)==1&&~isempty(cosx)==1&&isempty(sinx)~=1
            tw(1,boinx)=tw(1,boinx)*-1;
            attachterm=0;
        end
        % sin求导
        if ~isempty(sx)==1
            tw(2*qi+1,sx)=1;
            tw(2*qi,sx)=0;
        end
        % cos求导
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
            % 先剔除偏导等于0的项
            tw(:,setdiff(cop,scinx))=0;

            boinx=intersect(sinx,cosx);
            cx=setdiff(cosx,boinx);
            sx=setdiff(sinx,boinx);
            % 判断是否同一个变量同时含有sin和cos
            if ~isempty(boinx)==1&&~isempty(cosx)==1&&isempty(sinx)~=1
                tw(1,boinx)=tw(1,boinx)*-1;
                attachterm=0;
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
                % 先剔除列为零的,即偏导等于0的项；全给它们置零
                tw(:,setdiff(comple,qx))=0;
                tw(2*n+1+pn,qx)=0;
            end
        end
    else
        sprintf('Wrong! May be that pn is greater than n.')
    end
    attachterm=0;
elseif deg==2 % cal C(q,dq)
    if pn==0 % 只有revolute
        sinx=find(tw(2*qi,:));
        cosx=find(tw(2*qi+1,:));

        scinx=union(sinx,cosx);

        comple=1:size(tw,2);
        % 先剔除列为零的,即偏导等于0的项；全给它们置零
        tw(:,setdiff(comple,scinx))=0;

        % 符合函数求导
        boinx=intersect(sinx,cosx);
        cx=setdiff(cosx,boinx);
        sx=setdiff(sinx,boinx);
        % 判断是否同一个变量同时含有sin和cos
        if ~isempty(boinx)==1&&~isempty(cosx)==1&&isempty(sinx)~=1
            % [sin(q)cos(q)]' 对q求导等于cos^2(q)-sin^2(q)
            attachter=tw;
            attachter(2*qi:2*qi+1,boinx)=0;
            attachterm=attachter(:,boinx);
            
            tw(1,boinx)=tw(1,boinx)*-2;
            tw(2*qi,boinx)=2;
            tw(2*qi,boinx)=0;
            
            
        end
        % sin求导
        if ~isempty(sx)==1
            sx1=(tw(2*qi,sx)==1);
            tw(2*qi+1,sx1)=1;
            tw(2*qi,sx1)=0;
            
            sx2=(tw(2*qi,sx)==2);
            tw(2*qi+1,sx2)=1;
            tw(2*qi,sx2)=1;
            tw(1,sx2)=2*tw(1,sx2);
        end
        % cos求导
        if ~isempty(cx)==1
            cx1=(tw(2*qi+1,cx)==1);
            tw(2*qi,cx1)=1;
            tw(2*qi+1,cx1)=0;
            tw(1,cx1)=-1*tw(1,cx1);
            
            cx2=(tw(2*qi+1,cx)==2);
            tw(2*qi,cx2)=1;
            tw(2*qi+1,cx2)=1;
            tw(1,cx2)=-2*tw(1,cx2);
        end 
    elseif pn~=0&&pn<=n % revolute and prismatrical
        if qi~=pn
            sinx=find(tw(2*qi,:));
            cosx=find(tw(2*qi+1,:));

            scinx=union(sinx,cosx);
            cop=1:size(tw,2);
            % 先剔除偏导等于0的项
            tw(:,setdiff(cop,scinx))=0;

            boinx=intersect(sinx,cosx);
            cx=setdiff(cosx,boinx);
            sx=setdiff(sinx,boinx);
            % 判断是否同一个变量同时含有sin和cos
            if ~isempty(boinx)==1&&~isempty(cosx)==1&&isempty(sinx)~=1
                attachter=tw;
                attachter(2*qi:2*qi+1,boinx)=0;
                attachterm=attachter(:,boinx);
            
                tw(1,boinx)=tw(1,boinx)*-2;
                tw(2*qi,boinx)=2;
                tw(2*qi,boinx)=0;                
            end
            if isempty(sx)~=1
                sx1=(tw(2*qi,sx)==1);
                tw(2*qi+1,sx1)=1;
                tw(2*qi,sx1)=0;

                sx2=(tw(2*qi,sx)==2);
                tw(2*qi+1,sx2)=1;
                tw(2*qi,sx2)=1;
                tw(1,sx2)=2*tw(1,sx2);
            end
            if isempty(cx)~=1
                cx1=(tw(2*qi+1,cx)==1);
                tw(2*qi,cx1)=1;
                tw(2*qi+1,cx1)=0;
                tw(1,cx1)=-1*tw(1,cx1);

                cx2=(tw(2*qi+1,cx)==2);
                tw(2*qi,cx2)=1;
                tw(2*qi+1,cx2)=1;
                tw(1,cx2)=-2*tw(1,cx2);
            end
        else
            qx=find(tw(2*n+1+pn,:));
            if isempty(qx)~=1
                comple=1:size(tw,2);
                % 先剔除列为零的,即偏导等于0的项；全给它们置零
                tw(:,setdiff(comple,qx))=0;
                tw(2*n+1+pn,qx)=0;
            end
        end
    else
        sprintf('Wrong! May be that pn is greater than n.')
    end
    
else
    sprintf('Wrong! May be that pn is greater than n.')
end
%     
