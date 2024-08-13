function Dij = cal_dij_num(i,j,n,traceV,dtraceW,rn)
%cal_Dij 
% Dmat is inertia matrix of robot manipulator, its computation is based on
% the formula in the book "机器人动力学与控制-霍伟.pdf"  
%
%
% 中文：None
%
% Input:
% 	i:the first loop index
%   j:the second loop index, along with i,which can get the index k
%   Tmat:transform matrix
%   n:number of joint angles
%   q:joint angle vector
%   traceV:traceV and traceW is calculated sperately for high efficiency.
%
%
% Output:
% Dmat: indirect D matrix, that is V for  vector which is made of row-wise
% vectorized D.
% totalV: it is about dynamic parameters, a vector.
% 
% Example:
% 	  None
%   
%


% 利用公式计算Dij,返回dij指针
% [q,~,~,~]=fetch_variable(n);
kk=max([i,j]);
dij_v=cell(length(kk:n),1);
dij_w=cell(length(kk:n),1);
% t1=0;
for t=kk:n

    [rr,trace_k_w]=cal_Trace(dtraceW,t,i,j);

    trace_k_v=double(rr).*traceV{t};%换符号方向
    dij_w{t}=trace_k_w;%
    dij_v{t}=trace_k_v;%lie向量堆叠
end
dij_w1=horzcat(dij_w{:});
dij_v1=horzcat(dij_v{:});
sprintf('dij!')

[kw,kv]=cancelRZ_num_best(dij_v1,dij_w1);
if isempty(kw)==1
    kw=zeros(3*n,1);
    kv=zeros(size(traceV{1},1),1);
end

if isempty(kw)==1
    Dij.W=zeros(3*n,1);
    Dij.V=zeros(size(traceV{1},1),1).';
else
    % 输出行为function个数，列为参数
    [Dij.W,Dij.V] = sin2cos2Sim_num(kv,kw,n);
%      Dij.W=kw;
%      Dij.V=kv.';
end
% ckecks(Dij.W,Dij.V,t1,n,q);

end

% local function
function [rs,kroninx] = cal_Trace(traceW,t,i,j)
%   求解单个齐次变换矩阵的trace
% twj=calderivite(traceW.tw,j,rn);
% twi=calderivite(traceW.tw,i,rn);
% 
% %这里应该将T编程一个sparse的kronecker积,相加时再扩张数组.
twj=traceW{j,t}.difftw;
twi=traceW{i,t}.difftw;
kroninx=encode_kron(twj(2:end,:),twi(2:end,:));
rs=kron(twj(1,:),twi(1,:));%行向量
end

