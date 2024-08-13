x_end= zeros(3,1);
x_start=[-1,-1.8,-1.9]';
v_start = zeros(3,1);
v_end = zeros(3,1);
v_max = 0.8;
a = 0.5;
d = 0.5;

for i =1:3
    qr = generate_trajectory_points(v_start(i,1), v_end(i,1), v_max, a, d, x_start(i,1), x_end(i,1));
%     qr_Trj(:,i)=qr(:,1);
end