function [q,dq,dqq,gsym] = fetch_variable(n)

syms q1 q2 q3 q4 q5 q6 q7 real;
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
syms dqq1 dqq2 dqq3 dqq4 dqq5 dqq6 dqq7 real;
syms g;
q=[q1,q2,q3,q4,q5,q6,q7];
dq=[dq1,dq2,dq3,dq4,dq5,dq6,dq7];
ddq=[dqq1,dqq2,dqq3,dqq4,dqq5,dqq6,dqq7];
q=q(1:n);
dq=dq(1:n);
dqq=ddq(1:n);
gsym=g;
end

