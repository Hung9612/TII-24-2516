function [dtw,ddtw] = cal_traceWnew(Tmat,q,pn)
%CALTRACEW 
n=length(q);
base=funcBase(n);
dtw=cell(n,n);
ddtw=cell(n,n,n);
for p=1:n
    %
    [rj,pjw]=conformSgn(Tmat{1,p}.W);
    %
    tw=encoders(pjw,base,rj,pn);
    for j=1:p
        dtw{j,p}.difftw=cal_derivite_best(tw,j,pn);%
        for k=1:j
            ddtw{k,j,p}.difftw=cal_derivite_best(dtw{j,p}.difftw,k,pn);
        end
    end
end
    sprintf('º∆À„traceW_difftraceW£°')

end

