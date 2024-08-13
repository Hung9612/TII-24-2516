function flag = JointLimit(q,dq,qUp,dqUP)
% 电机最大转速为2000r/min；
% 关节角速度= 2000/160/60*2*pi =1.3090
% rad, rad/s
if abs(q(1))<qUp && abs(q(2))<qUp && abs(q(3))<qUp ...
    && abs(dq(1))<dqUP && abs(dq(2))<dqUP && abs(dq(3))<dqUP
    flag = 1;
else 
    flag = 0;
end
end

