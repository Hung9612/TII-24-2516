function Wmatrix = W_tau(in1,in2,in3)
%W_TAU
%    WMATRIX = W_TAU(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    13-Aug-2024 21:56:30

ddq1 = in3(1,:);
ddq2 = in3(2,:);
ddq3 = in3(3,:);
dq1 = in2(1,:);
dq2 = in2(2,:);
dq3 = in2(3,:);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = sin(q1);
t6 = sin(q2);
t7 = sin(q3);
t8 = dq1.^2;
t9 = dq2.^2;
t10 = dq3.^2;
t11 = q2.*2.0;
t12 = q3.*2.0;
t26 = ddq1.*6.5536e+4;
t13 = cos(t11);
t14 = t3.^2;
t15 = t4.^2;
t16 = sin(t11);
t17 = sin(t12);
t18 = t7.^2;
t19 = ddq1.*t3;
t20 = ddq2.*t7;
t22 = ddq1.*t4.*t6;
t24 = t4.*t8;
t27 = dq1.*dq2.*t3.*t6.*2.0;
t30 = t5.*6.59393e+5;
t31 = t3.*t6.*t7.*t8;
t21 = t7.*t19;
t23 = ddq1.*t14;
t25 = -t22;
t28 = -t24;
t29 = -t27;
t32 = (t8.*t16)./2.0;
t33 = t26+t30;
Wmatrix = reshape([ddq2.*t3-t6.*t9,t19,0.0,ddq1.*t4.*-2.0-t4.*t5.*1.006153869628906e+1+t6.*t20+dq1.*dq3.*t7.*2.0-ddq3.*t3.*t4-t2.*t3.*t7.*1.006153869628906e+1+t3.*t7.*t9+t3.*t7.*t10+dq2.*dq3.*t4.*t6.*2.0,(t6.*t7.*t33)./6.5536e+4,t2.*t7.*(-1.006153869628906e+1)-t7.*t8-t4.*t19-t3.*t4.*t5.*1.006153869628906e+1,t27+ddq1.*t6.^2,ddq2-t32,0.0,t3.*t20+dq1.*dq2.*t4.*2.0-ddq3.*t4.*t6-t6.*t7.*t9+t6.*t7.*t10-t4.*t6.*t19.*2.0-dq1.*dq2.*t4.*t14.*4.0+dq1.*dq3.*t3.*t6.*t7.*2.0,t21+t28+ddq3.*t7+t4.*t10+t14.*t24.*2.0+dq1.*dq3.*t3.*t4.*2.0,t20+t25-t31-dq1.*dq2.*t3.*t4.*2.0,t6.*t21-dq1.*dq2.*t7+(ddq2.*t3.*t4)./2.0+(ddq3.*t6.*t7)./2.0-(t4.*t6.*t9)./2.0+(t4.*t6.*t10)./2.0+dq1.*dq2.*t7.*t14.*2.0+dq1.*dq3.*t3.*t4.*t6,(ddq3.*t4)./2.0+(t7.*t8)./2.0-(t7.*t10)./2.0+(t4.*t19)./2.0-t7.*t8.*t14-dq1.*dq3.*t3.*t7,(ddq2.*t4)./2.0+(ddq1.*t6.*t7)./2.0-(t3.*t6.*t24)./2.0+dq1.*dq2.*t3.*t7,-ddq3.*t3-ddq1.*t15-t18.*t23+t18.*t27+dq1.*dq3.*t17+t4.*t6.*t20+dq2.*dq3.*t6.*t15.*2.0+t3.*t4.*t7.*t9-dq1.*dq3.*t4.*t7.*t14.*2.0,-t7.*(t20+t25+t31+dq2.*dq3.*t4.*2.0+dq1.*dq3.*t6.*t7.*2.0),-ddq3-t19-(t8.*t17)./2.0+(t9.*t17)./2.0+t7.*t14.*t24+dq1.*dq2.*t6.*t18.*2.0,dq1.*dq3+(ddq2.*t6)./2.0-(ddq1.*t17)./2.0+(t3.*t9)./2.0-dq1.*dq3.*t14-dq1.*dq3.*t15.*2.0-ddq2.*t6.*t15-t3.*t9.*t15+t4.*t7.*t23+dq1.*dq3.*t14.*t15.*2.0+dq2.*dq3.*t4.*t6.*t7.*2.0-dq1.*dq2.*t3.*t4.*t6.*t7.*2.0,-dq2.*dq3+(ddq1.*t6)./2.0+t4.*t20+dq2.*dq3.*t15.*2.0-ddq1.*t6.*t15+t3.*t6.*t7.*t24+dq1.*dq3.*t4.*t6.*t7.*2.0,t8.*(-1.0./2.0)+t9./2.0+(t8.*t14)./2.0+t8.*t15-t9.*t15-t8.*t14.*t15-dq1.*dq2.*t4.*t6.*t7.*2.0,ddq1+t23+t29+ddq3.*t3.*2.0-dq2.*dq3.*t6.*2.0,ddq2+t32+dq1.*dq3.*t6.*2.0,ddq3.*2.0+t19.*2.0-dq1.*dq2.*t6.*2.0,(ddq1.*t16)./2.0+dq1.*dq2.*t13,t8.*t13.*(-1.0./2.0),0.0,t23+t29,ddq2+t32,0.0,ddq1.*t7.*2.0+t5.*t7.*1.006153869628906e+1+dq1.*dq3.*t4.*2.0+ddq2.*t4.*t6+ddq3.*t3.*t7-t2.*t3.*t4.*1.006153869628906e+1+t3.*t4.*t9+t3.*t4.*t10-dq2.*dq3.*t6.*t7.*2.0,(t4.*t6.*t33)./6.5536e+4,t21+t28-t2.*t4.*1.006153869628906e+1+t3.*t5.*t7.*1.006153869628906e+1,ddq2.*t6+t3.*t9,ddq1.*t6,0.0,ddq1,0.0,0.0,t2,0.0,0.0,t2.*t6,t3.*t5,0.0,t2.*t3,-t5.*t6,0.0,t5,0.0,0.0],[3,17]);
