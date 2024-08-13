function [dtw,ddtw] = cal_traceWnew(Tmat,q,pn)
%CALTRACEW 
n=length(q);
base=funcBase(n);
dtw=cell(n,n);
ddtw=cell(n,n,n);
for p=1:n
    % ����ÿ����α任����ĺ�������
    [rj,pjw]=conformSgn(Tmat{1,p}.W);
    % �Ե�����α任������б���
    tw=encoders(pjw,base,rj,pn);
    for j=1:p
        dtw{j,p}.difftw=cal_derivite_best(tw,j,pn);% �Ե�p���任���������еĵ���
        for k=1:j
            ddtw{k,j,p}.difftw=cal_derivite_best(dtw{j,p}.difftw,k,pn);
        end
    end
end
    sprintf('����traceW_difftraceW��')

end

