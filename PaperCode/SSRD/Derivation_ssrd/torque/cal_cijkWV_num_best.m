function cijk = cal_cijkWV_num_best(dtraceW,ddtracew,traceV,i,j,k,n)

kk=max([i,j,k]);
cijk_v=cell(length(kk:n),1);
cijk_w=cell(length(kk:n),1);
for t=kk:n

    [rr,trace_k_w]=cal_Trace(dtraceW,ddtracew,t,i,j,k);
    trace_k_v=double(rr).*traceV{t};
    cijk_w{t}=trace_k_w;
    cijk_v{t}=trace_k_v;
end
cijk_w1=horzcat(cijk_w{:});
cijk_v1=horzcat(cijk_v{:});
% sprintf('cijk!')
[kw,kv]=cancelRZ_num_best(cijk_v1,cijk_w1);

if isempty(kw)==1
    cijk.W=zeros(3*n,1);
    cijk.V=zeros(size(traceV{1},1),1).';
else
    
    [cijk.W,cijk.V] = sin2cos2Sim_num(kv,kw,n);
end

end

% local function
function [rs,kroninx] = cal_Trace(dtraceW,ddtraceW,t,i,j,k)

twkj=ddtraceW{j,k,t}.difftw;
twi=dtraceW{i,t}.difftw;
kroninx=encode_kron(twi(2:end,:),twkj(2:end,:));
rs=kron(twi(1,:),twkj(1,:));
end
