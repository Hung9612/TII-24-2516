function [tauw,P,Theta,nulltheta] = cal_base_num(taud,tauc,taug,opt)

if opt==2
    tol=9;
    X=cal_Yd(opt);
    Jmat=[taud.v;tauc.v;taug.v];
    [Lp,Up]=FRdecomposition_SVD_thin(Jmat,opt);

    Z=[Up(:,1:20);
       Up(:,21:40)];
    [Lz,Uz]=FRdecomposition_SVD_thin(Z,opt);

    Theta=vpa(Uz*X.',9);
    nulltheta=X(any(round(Uz,tol)));
    sprintf('The number of base parameter is %d',length(Theta))

    %
    ro=size(Up,1);
    P1=Lp*Lz(1:ro,:);
    P2=Lp*Lz(ro+1:2*ro,:);
   
    P={round(P1,tol),round(P2,tol)};
    tauw=double([taud.w,tauc.w,taug.w]);
    
    elseif opt==3
    tol=9;
    X=cal_Yd(opt);
    Jmat=[taud.v;tauc.v;taug.v];
    [Lp,Up]=FRdecomposition_SVD_thin(single(Jmat),opt);
    Z=[Up(:,1:30);
       Up(:,31:60);
       Up(:,61:90)];
    [Lz,Uz]=FRdecomposition_SVD_thin(Z,opt);

    Theta=Uz*X.';
    nulltheta=X(any(round(Uz,tol)));
    sprintf('The number of base parameter is %d',length(Theta))

    %
    ro=size(Up,1);
    P1=Lp*Lz(1:ro,:);
    P2=Lp*Lz(ro+1:2*ro,:);
    P3=Lp*Lz(2*ro+1:3*ro,:);

    P={P1,P2,P3};

    tauw=double([taud.w,tauc.w,taug.w]);
    
   elseif opt==4
    tol=9;
    X=cal_Yd(opt);
    Jmat=[taud.v;tauc.v;taug.v];
    [Lp,Up]=FRdecomposition_SVD_thin(single(Jmat),opt);
    Z=[Up(:,1:40);
       Up(:,41:80);
       Up(:,81:120);
       Up(:,121:160)];
    [Lz,Uz]=FRdecomposition_SVD_thin(Z,opt);

    Theta=vpa(Uz*X.',9);
    nulltheta=X(any(round(Uz,tol)));
    sprintf('The number of base parameter is %d',length(Theta))

  
    ro=size(Up,1);
    P1=Lp*Lz(1:ro,:);
    P2=Lp*Lz(ro+1:2*ro,:);
    P3=Lp*Lz(2*ro+1:3*ro,:);
    P4=Lp*Lz(3*ro+1:4*ro,:);
   
    P={round(P1,tol),round(P2,tol),round(P3,tol),round(P4,tol)};
    tauw=double([taud.w,tauc.w,taug.w]);

elseif opt==7
    tol=9;
    X=cal_Yd(opt);
    Jmat=[taud.v;tauc.v;taug.v];

    [Lp,Up]=FRdecomposition_SVD_thin(single(Jmat),opt);

    Z=[Up(:,1:70);
       Up(:,71:140);
       Up(:,141:210);
       Up(:,211:280);
       Up(:,281:350);
       Up(:,351:420);
       Up(:,421:490)];

    [Lz,Uz]=FRdecomposition_SVD_thin(Z,opt);
    Theta=vpa(Uz*X.',9);%最小惯性参数集
    nulltheta=X(sum(abs(round(Uz,tol)),1)==0);
    sprintf('The number of base parameter is %d',length(Theta))

    %
    ro=size(Up,1);
    P1=Lp*Lz(1:ro,:);
    P2=Lp*Lz(ro+1:2*ro,:);
    P3=Lp*Lz(2*ro+1:3*ro,:);
    P4=Lp*Lz(3*ro+1:4*ro,:);
    P5=Lp*Lz(4*ro+1:5*ro,:);
    P6=Lp*Lz(5*ro+1:6*ro,:);
    P7=Lp*Lz(6*ro+1:7*ro,:);
    P={round(P1,tol),round(P2,tol),round(P3,tol),round(P4,tol),round(P5,tol),round(P6,tol),round(P7,tol)};
    tauw=double([taud.w,tauc.w,taug.w]);

    
elseif opt==5
    tol=9;
    X=cal_Yd(opt);
    Jmat=[taud.v;tauc.v;taug.v];
    [Lp,Up]=FRdecomposition_SVD_thin(Jmat,opt);

    Z=[Up(:,1:50);
       Up(:,51:100);
       Up(:,101:150);
       Up(:,151:200);
       Up(:,201:250)];

    [Lz,Uz]=FRdecomposition_SVD_thin(Z,opt);

    Theta=vpa(Uz*X.',9);
    nulltheta=X(sum(abs(round(Uz,tol)),1)==0);
    sprintf('The number of base parameter is %d',length(Theta))
    
    ro=size(Up,1);
    P1=Lp*Lz(1:ro,:);
    P2=Lp*Lz(ro+1:2*ro,:);
    P3=Lp*Lz(2*ro+1:3*ro,:);
    P4=Lp*Lz(3*ro+1:4*ro,:);
    P5=Lp*Lz(4*ro+1:5*ro,:);
    P={round(P1,tol),round(P2,tol),round(P3,tol),round(P4,tol),round(P5,tol)};
    tauw=double([taud.w,tauc.w,taug.w]);    
elseif opt==6
    tol=9;
    X=cal_Yd(opt);
    Jmat=[taud.v;tauc.v;taug.v];
    [Lp,Up]=FRdecomposition_SVD_thin(Jmat,opt);

    Z=[Up(:,1:60);
       Up(:,61:120);
       Up(:,121:180);
       Up(:,181:240);
       Up(:,241:300);
       Up(:,301:360)];

    [Lz,Uz]=FRdecomposition_SVD_thin(Z,opt);

    Theta=vpa(Uz*X.',9);
    nulltheta=X(sum(abs(round(Uz,tol)),1)==0);
    sprintf('The number of base parameter is %d',length(Theta))
    
    ro=size(Up,1);
    P1=Lp*Lz(1:ro,:);
    P2=Lp*Lz(ro+1:2*ro,:);
    P3=Lp*Lz(2*ro+1:3*ro,:);
    P4=Lp*Lz(3*ro+1:4*ro,:);
    P5=Lp*Lz(4*ro+1:5*ro,:);
    P6=Lp*Lz(5*ro+1:6*ro,:);
    P={round(P1,tol),round(P2,tol),round(P3,tol),round(P4,tol),round(P5,tol),round(P6,tol)};
    tauw=double([taud.w,tauc.w,taug.w]);

else
    sprintf('error of opt!')
end

end

