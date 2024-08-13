function tauC_i = Decode_TauC(tauC,t)
syms q1 q2 q3 dq1 dq2 dq3
w=tauC.w;
v=tauC.v;
NN = length(w(1,:));
phi = w;
w_bar=zeros(1,NN)*q1;
for i=1:NN
    w_bar(i) = (sin(q1))^phi(1,i)*(cos(q1))^phi(2,i)...
               *(sin(q2))^phi(3,i)*(cos(q2))^phi(4,i)...
               *(sin(q3))^phi(5,i)*(cos(q3))^phi(6,i)...
               *(dq1)^phi(10,i)*(dq2)^phi(11,i)*(dq3)^phi(12,i);
end
X = cal_Yd(3);
if t==1
   tauC_i = w_bar *v*kron([1;0;0],X.');
elseif t==2
   tauC_i = w_bar *v*kron([0;1;0],X.');
else
   tauC_i = w_bar *v*kron([0;0;1],X.');

end
end

