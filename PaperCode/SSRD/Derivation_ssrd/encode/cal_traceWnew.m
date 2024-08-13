function [dtw,ddtw] = cal_traceWnew(Tmat,q,pn)
%CALTRACEW 
n=length(q);
base=funcBase(n);
dtw=cell(n,n);
ddtw=cell(n,n,n);
for p=1:n
    % 计算每个齐次变换矩阵的函数向量
    [rj,pjw]=conformSgn(Tmat{1,p}.W);
    % 对单个齐次变换矩阵进行编码
    tw=encoders(pjw,base,rj,pn);
    for j=1:p
        dtw{j,p}.difftw=cal_derivite_best(tw,j,pn);% 对第p个变换矩阵求所有的导数
        for k=1:j
            ddtw{k,j,p}.difftw=cal_derivite_best(dtw{j,p}.difftw,k,pn);
        end
    end
end
    sprintf('计算traceW_difftraceW！')

end

