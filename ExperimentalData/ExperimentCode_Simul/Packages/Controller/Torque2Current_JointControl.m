function [Current,Torq_jointSpa] = Torque2Current_JointControl(u_PD,u_Ada,uff,inputSat, HomeFlag)
acc1 = 0.001*[7, 6.21, 2.39].';
B3 = [160, 160, 100].';
Ki = [0.194, 0.109, 0.111].';

%  Torq_jointSpa = Torq.*[1; 1; -1] ;

%  Torq_jointSpa = Torq.*[1; 1; -1] + uff;
if HomeFlag ==1
%      Torq_jointSpa = (u_PD ).*[1; 1; -1] + uff;
          Torq_jointSpa =  uff;

else
%      Torq_jointSpa = (u_PD.*[1; 1; -1] + u_Ada.*[1; 1; 1]) + uff;
     Torq_jointSpa = u_PD.*[1; 1; -1] + u_Ada.*[1; 1; 1]; 
%      Torq_jointSpa = u_PD.*[1; 1; -1]; 

end
if inputSat(1) < abs(Torq_jointSpa(1))
    Torq_jointSpa(1) = inputSat(1);
end

if inputSat(2) < abs(Torq_jointSpa(2))
    Torq_jointSpa(2) = inputSat(2);
end

if inputSat(3) < abs(Torq_jointSpa(3))
    Torq_jointSpa(3) = inputSat(3);
end

Current = Torq_jointSpa./Ki./B3./acc1;
end

