function [tw,dtw] = cal_tracewEncode(Tmat,pn)
%CALTRACEW 
%
n=length(Tmat);
base=funcBase(n);
% tw=cell(1,n);
dtw=cell(1,n);

for p=1:n
    %
        [rj,pjw]=conformSgn(Tmat{1,p}.W);
        %
%         tw{1,p}.tw=encoders(pjw,base,rj,pn);
        tw=encoders(pjw,base,rj,pn);
    for j=1:p
        dtw{j,p}.difftw=cal_derivite(tw,j,pn);%
    end
end

end

