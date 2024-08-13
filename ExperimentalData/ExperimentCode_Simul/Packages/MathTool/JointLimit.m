function flag = JointLimit(q,dq,qUp,dqUP)
% 
% 
% rad, rad/s
if abs(q(1))<qUp && abs(q(2))<qUp && abs(q(3))<qUp ...
    && abs(dq(1))<dqUP && abs(dq(2))<dqUP && abs(dq(3))<dqUP
    flag = 1;
else 
    flag = 0;
end
end

