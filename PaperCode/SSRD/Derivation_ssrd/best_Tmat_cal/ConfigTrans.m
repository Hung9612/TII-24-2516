function [T,E] = ConfigTrans(theta)
syms d1;syms d2;syms d3; syms d4;
syms z1;syms z2;syms z3; syms z4; syms z5; syms z6;syms z7;

n=length(theta);
if n==6
        q1=theta(1);q2=theta(2);q3=theta(3);q4=theta(4);q5=theta(5);q6=theta(6);

    e1=[zrot(q1),zeros(1,3)';0,0,0,1];
    e2=[yrot(q2),zeros(1,3)';0,0,0,1];
    e3=[yrot(q3),zeros(1,3)';0,0,0,1];
    e4=[yrot(q4),zeros(1,3)';0,0,0,1];
    e5=[zrot(q5),zeros(1,3)';0,0,0,1];
    e6=[yrot(q6),zeros(1,3)';0,0,0,1];
    T6=transl([0,0,0.08]);
    T5=transl([0,0,0.365]);
    T4=transl([0,0,0.035]);
    T3=transl([0,0,0.315]);
    T2=transl([0.025,0,0.183]);
    T1=transl([0,0,0]);
    
    T={T1,T2,T3,T4,T5,T6};
    E = {e1,e2,e3,e4,e5,e6};
elseif n==7
% sym type
   q1=theta(1);q2=theta(2);
   q3=theta(3);q4=theta(4);
   q5=theta(5);q6=theta(6);q7=theta(7);
    
    e1=[zrot(q1),zeros(1,3)';0,0,0,1];
    e2=[xrot(q2),zeros(1,3)';0,0,0,1];
    e3=[zrot(q3),zeros(1,3)';0,0,0,1];
    e4=[xrot(q4),zeros(1,3)';0,0,0,1];
    e5=[zrot(q5),zeros(1,3)';0,0,0,1];
    e6=[xrot(q6),zeros(1,3)';0,0,0,1];
    e7=[zrot(q7),zeros(1,3)';0,0,0,1];

    T7=transl([0,0,0.165]);
    T6=transl([0,0,0.165]);
    T5=transl([0,0,0.185]);
    T4=transl([0,0,0.165]);
    T3=transl([0,0,0.185]);
    T2=transl([0,0,0.180]);
    T1=transl([0,0,0.1435]);
   T={T1,T2,T3,T4,T5,T6,T7};
    E = {e1,e2,e3,e4,e5,e6,e7};

    
elseif n==4%SCARA
   q1=theta(1);q2=theta(2);
   q3=theta(3);q4=theta(4);
   
    
    e1=[zrot(q1),zeros(1,3)';0,0,0,1];
    e2=[zrot(q2),zeros(1,3)';0,0,0,1];
    e3=[eye(3),[0,0,-q3].';0,0,0,1];
    e4=[zrot(q4),zeros(1,3)';0,0,0,1];

    T4=transl([0,0,0]);
    T3=transl([0.26,0,0]);
    T2=transl([0.45,0,0]);
    T1=transl([0,0,0]);
   T={T1,T2,T3,T4};
    E = {e1,e2,e3,e4};
elseif n==2
   q1=theta(1);q2=theta(2);

    
    e1=[zrot(q1),zeros(1,3)';0,0,0,1];
    e2=[zrot(q2),zeros(1,3)';0,0,0,1];

    T2=transl([0.25,0,0]);
    T1=transl([0,0,0]);
   T={T1,T2};
    E = {e1,e2};
elseif n==5
% sym type
   q1=theta(1);q2=theta(2);
   q3=theta(3);q4=theta(4);
   q5=theta(5);
    
    e1=[zrot(q1),zeros(1,3)';0,0,0,1];
    e2=[xrot(q2),zeros(1,3)';0,0,0,1];
    e3=[zrot(q3),zeros(1,3)';0,0,0,1];
    e4=[xrot(q4),zeros(1,3)';0,0,0,1];
    e5=[zrot(q5),zeros(1,3)';0,0,0,1];
 


    T5=transl([0,0,0.185]);
    T4=transl([0,0,0.165]);
    T3=transl([0,0,0.185]);
    T2=transl([0,0,0.180]);
    T1=transl([0,0,0.1435]);
   T={T1,T2,T3,T4,T5};
    E = {e1,e2,e3,e4,e5};
    
%     elseif n==3
% % sym type
%    q1=theta(1);q2=theta(2);
%    q3=theta(3);
%     
%     e1=[zrot(q1),zeros(1,3)';0,0,0,1];
%     e2=[xrot(q2),zeros(1,3)';0,0,0,1];
%     e3=[zrot(q3),zeros(1,3)';0,0,0,1];
%     T3=transl([0,0,0.185]);
%     T2=transl([0,0,0.180]);
%     T1=transl([0,0,0.1435]);
%    T={T1,T2,T3};
%     E = {e1,e2,e3};
    elseif n==3
% sym type
   q1=theta(1);q2=theta(2);
   q3=theta(3);
    
    e1=[xrot(q1),zeros(1,3)';0,0,0,1];
    e2=[zrot(q2),zeros(1,3)';0,0,0,1];
    e3=[xrot(q3),zeros(1,3)';0,0,0,1];
    T0 = transl([0,0,0.3245]);
    T1 = transl([0,0,0.35]);
    T2 = transl([0,0,0.35]);
    T3 = transl([0,0,0.275]);
    T = {T0,T1,T2*T3};
    E = {e1,e2,e3};

else
    sprintf('enable to calculate 4,6,7 joints! Please input other angles.');
end
end
