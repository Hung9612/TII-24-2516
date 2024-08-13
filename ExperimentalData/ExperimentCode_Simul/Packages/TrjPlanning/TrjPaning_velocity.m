function [qr1,qr2,qr3] = TrjPaning_velocity(Q_end,Q_start)
v_start = 0;
v_end = 0;
v_max = 0.8;
a = 0.5;
d = 0.5;
for i=1:3
flag = 0;
    if Q_end(i)<Q_start(i)
        StartP = Q_end(i);
        EndPf  = Q_start(i);
        flag=1;
    else
        StartP = Q_start(i);
        EndPf  = Q_end(i);
    end
    
qr_i = generate_trajectory_points(v_start, v_end, v_max, a, d, StartP, EndPf);
if flag == 1
    Qr{i} = flip(qr_i);
else
    Qr{i} = qr_i;
end
% for i =1:3
%     qr = generate_trajectory_points(v_start(i,1), v_end(i,1), v_max, a, d, x_start(i,1), x_end(i,1));
% %     qr_Trj(:,i)=qr(:,1);
end
qr1 = Qr{1};
qr2 = Qr{2};
qr3 = Qr{3};


