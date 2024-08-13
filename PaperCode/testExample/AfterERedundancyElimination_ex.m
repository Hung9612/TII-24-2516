
rbt_n.name='DOF3_Ctrl';
rbt_n.DoF=3;
rbt_n.pn=0;
rbt_n.opt=1;
rbt_n.Gdir=[0,0,-1];
%%

[q,dq,dqq,~]=fetch_variable(rbt_n.DoF);
rbt_n.Tmat = cal_Tmat_best(q);


%%

[rbt_n.dmat,dtraceW,ddtracew,traceV]=cal_dmat(q,rbt_n.Tmat,rbt_n.opt,rbt_n.pn);


%%
rbt_n.taud = taucd_num(rbt_n.dmat,dqq);

%%
n=3;
for i=1:n
    error = simplify(TauD(i) - Decode_TauD(rbt_n.taud,i))
end

%%
[~,rbt_n.tauc] = taucc_num_best(dq,dtraceW,ddtracew,traceV,rbt_n.pn);
%%
n=3;
for i=1:n
    error = vpa(simplify(TauC(i) - Decode_TauC(rbt_n.tauc,i)),6)
end

%%
[rbt_n.taug]= taucg_num_fast(rbt_n.Tmat,dtraceW,q,rbt_n.Gdir);
%%
n=3;
for i=1:n
    error = simplify(TauG(i) - Decode_TauG(rbt_n.taug,i))
end


%%
resultRobot_n = reshape2(rbt_n);

%%
[rbt_n.tauw,rbt_n.Pmap,rbt_n.Theta,rbt_n.NullTheta] = cal_base_num(resultRobot_n.taud,resultRobot_n.tauc,resultRobot_n.taug,rbt_n.DoF);

%%
syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3
Q = [q1 q2 q3].';
dQ = [dq1 dq2 dq3].';
ddQ = [ddq1 ddq2 ddq3].';
NN = length(rbt_n.tauw(1,:));
phi = rbt_n.tauw;
w_bar=zeros(1,NN)*q1;
for i=1:NN
    w_bar(i) = (sin(q1))^phi(4,i)*(sin(q2))^phi(5,i)*(sin(q3))^phi(6,i)*(cos(q1))^phi(7,i)*(cos(q2))^phi(8,i)*(cos(q3))^phi(9,i)...
                *dq1^phi(10,i)*dq2^phi(11,i)*dq3^phi(12,i)*ddq1^phi(13,i)*ddq2^phi(14,i)*ddq3^phi(15,i);
end

for i=1:n
    error = vpa(simplify(Tau(i) - w_bar * rbt_n.Pmap{i}*rbt_n.Theta), 6);
    mapSymType(error, 'vpareal', @(x) piecewise(abs(x)<=1e-6, 0, x))
end

%%
Wmatrix = zeros(rbt_n.DoF,length(rbt_n.Theta))*q1;
for i=1:rbt_n.DoF
    Wmatrix(i,:) = simplify(w_bar*rbt_n.Pmap{1,i});
end
error = simplify(Tau - Wmatrix*rbt_n.Theta);
error = vpa(error,6); %
mapSymType(error, 'vpareal', @(x) piecewise(abs(x)<=1e-6, 0, x))
Re_GenW(Wmatrix,Q,dQ,ddQ)
Ccode_generte(rbt_n);

