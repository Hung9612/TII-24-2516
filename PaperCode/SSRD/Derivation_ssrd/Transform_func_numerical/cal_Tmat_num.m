function T = cal_Tmat_num(theta,numericals)

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
    
    T56=T6*e6;
    T45=T5*e5;
    T34=T4*e4;
    T23=T3*e3;
    T12=T2*e2;
    
    t01=T1*e1;
    t02=t01*T12;
    t03=t02*T23;
    t04=t03*T34;
    t05=t04*T45;
    t06=t05*T56;
    % return 
    % eg:T06.V_fast,T06.W
    [~,T06]=MatMultiply_best_num({T1,T2,T3,T4,T5,T6},{e1,e2,e3,e4,e5,e6},2);
    [~,T05]=MatMultiply_best_num({T1,T2,T3,T4,T5},{e1,e2,e3,e4,e5},2);
    [~,T04]=MatMultiply_best_num({T1,T2,T3,T4},{e1,e2,e3,e4},2);
    [~,T03]=MatMultiply_best_num({T1,T2,T3},{e1,e2,e3},2);
    [~,T02]=MatMultiply_num({T1,T2},{e1,e2},2);
    [~,T01]=MatMultiply_num({T1},{e1},2);
    
    erorT6=reshape(T06.V_fast*T06.W,[4,4]).'-t06;
    erorT5=reshape(T05.V_fast*T05.W,[4,4]).'-t05;
    erorT4=reshape(T04.V_fast*T04.W,[4,4]).'-t04;
    erorT3=reshape(T03.V_fast*T03.W,[4,4]).'-t03;
    erorT2=reshape(T02.V_fast*T02.W,[4,4]).'-t02;
    erorT1=reshape(T01.V_fast*T01.W,[4,4]).'-t01;
    
    error6=simplify(erorT6);
    error5=simplify(erorT5);
    error4=simplify(erorT4);
    error3=simplify(erorT3);
    error2=simplify(erorT2);
    error1=simplify(erorT1);
    if (error1+error2+error3+error3+error4+error5+error6)==0
        T={T01,T02,T03,T04,T05,T06};
    end
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
    % Jmat=subs(Jmat,{z1,z2,z3,z4,z5,z6,z7},{0.1435,0.180,0.185,0.165,0.185,0.165,0.165});

    T7=transl([0,0,0.165]);
    T6=transl([0,0,0.165]);
    T5=transl([0,0,0.185]);
    T4=transl([0,0,0.165]);
    T3=transl([0,0,0.185]);
    T2=transl([0,0,0.180]);
    T1=transl([0,0,0.1435]);
    
    T67=T7*e7;
    T56=T6*e6;
    T45=T5*e5;
    T34=T4*e4;
    T23=T3*e3;
    T12=T2*e2;
    t01=T1*e1;
    
    t02=t01*T12;
    t03=t02*T23;
    t04=t03*T34;
    t05=t04*T45;
    t06=t05*T56;
    t07=t06*T67;
    % return 
    % eg:T06.V_fast,T06.W
    [~,T07]=MatMultiply_best_num({T1,T2,T3,T4,T5,T6,T7},{e1,e2,e3,e4,e5,e6,e7},2);
    [~,T06]=MatMultiply_best_num({T1,T2,T3,T4,T5,T6},{e1,e2,e3,e4,e5,e6},2);
    [~,T05]=MatMultiply_best_num({T1,T2,T3,T4,T5},{e1,e2,e3,e4,e5},2);
    [~,T04]=MatMultiply_best_num({T1,T2,T3,T4},{e1,e2,e3,e4},2);
%     [~,T03]=MatMultiply_best_num({T1,T2,T3},{e1,e2,e3},2);
    [~,T03]=MatMultiply_num({T1,T2,T3},{e1,e2,e3},2);

    [~,T02]=MatMultiply_num({T1,T2},{e1,e2},2);
    [~,T01]=MatMultiply_num({T1},{e1},2);
    
    erorT7=reshape(T07.V_fast*T07.W,[4,4]).'-t07;
    erorT6=reshape(T06.V_fast*T06.W,[4,4]).'-t06;
    erorT5=reshape(T05.V_fast*T05.W,[4,4]).'-t05;
    erorT4=reshape(T04.V_fast*T04.W,[4,4]).'-t04;
    erorT3=reshape(T03.V_fast*T03.W,[4,4]).'-t03;
    erorT2=reshape(T02.V_fast*T02.W,[4,4]).'-t02;
    erorT1=reshape(T01.V_fast*T01.W,[4,4]).'-t01;
    
    error7=simplify(erorT7);
    error6=simplify(erorT6);
    error5=simplify(erorT5);
    error4=simplify(erorT4);
    error3=simplify(erorT3);
    error2=simplify(erorT2);
    error1=simplify(erorT1);
    if (error1+error2+error3+error4+error5+error6+error7)==0
        T={T01,T02,T03,T04,T05,T06,T07};
    end
    
elseif n==4%SCARA
    q1=theta(1);
    q2=theta(2);
    q3=theta(3);
    q4=theta(4);
    e1=[zrot(q1),zeros(1,3)';0,0,0,1];
    e2=[zrot(q2),zeros(1,3)';0,0,0,1];
    e3=[eye(3),[0,0,-q3].';0,0,0,1];
    e4=[zrot(q4),zeros(1,3)';0,0,0,1];
%     T4=transl([0,0,0]);
%     T3=transl([z2,0,0]);
%     T2=transl([z1,0,0]);
%     T1=transl([0,0,0]);
    T4=transl([0,0,0]);
    T3=transl([0.26,0,0]);
    T2=transl([0.45,0,0]);
    T1=transl([0,0,0]);
    T34=T4*e4;
    T23=T3*e3;
    T12=T2*e2;
    t01=T1*e1;
    t02=t01*T12;
    t03=t02*T23;
    t04=t03*T34;
    % return 
    % eg:T04.V_fast,T04.W
    [~,T04]=MatMultiply_num({T1,T2,T3,T4},{e1,e2,e3,e4},2);
    [~,T03]=MatMultiply_num({T1,T2,T3},{e1,e2,e3},2);
    [~,T02]=MatMultiply_num({T1,T2},{e1,e2},2);
    [~,T01]=MatMultiply_num({T1},{e1},2);

    erorT4=reshape(T04.V_fast*T04.W,[4,4]).'-t04;
    erorT3=reshape(T03.V_fast*T03.W,[4,4]).'-t03;
    erorT2=reshape(T02.V_fast*T02.W,[4,4]).'-t02;
    erorT1=reshape(T01.V_fast*T01.W,[4,4]).'-t01;

    error1=simplify(erorT4);
    error2=simplify(erorT3);
    error3=simplify(erorT2);
    error4=simplify(erorT1);
    if (error1+error2+error3+error3+error4)==0
        T={T01,T02,T03,T04};
    end

elseif n==2
    q1=theta(1);
    q2=theta(2);
    e1=[zrot(q1),zeros(1,3)';0,0,0,1];
    e2=[zrot(q2),zeros(1,3)';0,0,0,1];

    T2=transl([0.25,0,0]);
    T1=transl([0,0,0]);
    T12=T2*e2;
    t01=T1*e1;
    t02=t01*T12;

    % return 
    % eg:T02.V_fast,T02.W
    [~,T02]=MatMultiply_num({T1,T2},{e1,e2},2);
    [~,T01]=MatMultiply_num({T1},{e1},2);


    erorT2=reshape(T02.V_fast*T02.W,[4,4]).'-t02;
    erorT1=reshape(T01.V_fast*T01.W,[4,4]).'-t01;

    error1=simplify(erorT1);
    error2=simplify(erorT2);
   
    if (error1+error2)==0
        T={T01,T02};
    end
    
else
    sprintf('enable to calculate 4,6,7 joints! Please input other angles.');
end

end

