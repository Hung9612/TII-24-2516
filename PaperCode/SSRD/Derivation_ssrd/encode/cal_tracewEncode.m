function [tw,dtw] = cal_tracewEncode(Tmat,pn)
%CALTRACEW 
%   编码齐次变换矩阵
n=length(Tmat);
base=funcBase(n);
% tw=cell(1,n);
dtw=cell(1,n);

for p=1:n
    % 计算每个齐次变换矩阵的函数向量
        [rj,pjw]=conformSgn(Tmat{1,p}.W);
        % 对单个齐次变换矩阵进行编码
%         tw{1,p}.tw=encoders(pjw,base,rj,pn);
        tw=encoders(pjw,base,rj,pn);
    for j=1:p
        dtw{j,p}.difftw=cal_derivite(tw,j,pn);% 对第p个变换矩阵求所有的导数
    end
end

end

