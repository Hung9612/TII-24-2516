function dij = Decode_D(dmat,t)
syms q1 q2 q3 
w=dmat{t,1}.W;
v=dmat{t,1}.V;
NN = length(w(1,:));
phi = w;
w_bar=zeros(1,NN)*q1;
for i=1:NN
    w_bar(i) =  (sin(q1))^phi(1,i)*(cos(q1))^phi(2,i)...
               *(sin(q2))^phi(3,i)*(cos(q2))^phi(4,i)...
               *(sin(q3))^phi(5,i)*(cos(q3))^phi(6,i);
end
X = cal_Yd(3);
dij = w_bar *v*X.';
end

