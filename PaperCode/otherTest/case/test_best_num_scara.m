clear
rbt_n.name='scarav';
rbt_n.DoF=4;
rbt_n.pn=3;
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
[~,rbt_n.tauc] = taucc_num_best(dq,dtraceW,ddtracew,traceV,rbt_n.pn);

%%
 rbt_n.taug = taucg_num_fast(rbt_n.Tmat,dtraceW,q,rbt_n.Gdir);
 
%%
resultRobot_n = reshape2(rbt_n);

%%
[rbt_n.tauw,rbt_n.Pmap,rbt_n.Theta,rbt_n.NullTheta] = ...
    cal_base_num(resultRobot_n.taud,resultRobot_n.tauc,...
    resultRobot_n.taug,rbt_n.DoF);Ccode_generte(rbt_n);
%% 
Ccode_generte(rbt_n);
%%
filenames = sprintf('cacheData_%s.mat',rbt_n.name);
save(filenames,'rbt_n');



