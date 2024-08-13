function out= funcBase(n)
%UNTITLED 
syms q1;syms q2;syms q3;syms q4;syms q5;syms q6;syms q7;

if n==2
    out=[sin(q1) ; cos(q1) ; sin(q2) ; cos(q2);q1; q2  ];
    elseif n==3
    out=[sin(q1) ; cos(q1); sin(q2) ; cos(q2) ;   ...
         sin(q3) ; cos(q3);q1; q2; q3  ];
elseif n==4
    out=[sin(q1) ; cos(q1) ;sin(q2) ; cos(q2) ; sin(q3) ; cos(q3) ;  ...
         sin(q4) ; cos(q4);q1; q2; q3; q4  ];
elseif n==7
    out=[sin(q1) ; cos(q1) ; sin(q2) ; cos(q2) ; sin(q3) ; cos(q3) ;...
         sin(q4) ; cos(q4) ; sin(q5) ; cos(q5) ; sin(q6) ; cos(q6) ;...
         sin(q7) ; cos(q7);q1; q2; q3; q4;q5; q6; q7 ];
elseif n==6
    out=[ sin(q1) ; cos(q1) ; sin(q2) ; cos(q2) ; sin(q3) ; cos(q3) ;...
         sin(q4) ; cos(q4) ; sin(q5) ; cos(q5) ; sin(q6) ; cos(q6);q1;...
         q2; q3; q4;q5; q6];
     elseif n==5
    out=[ sin(q1) ; cos(q1) ; sin(q2) ; cos(q2) ; sin(q3) ; cos(q3) ;...
         sin(q4) ; cos(q4) ; sin(q5) ; cos(q5) ;  q1;...
         q2; q3; q4;q5];
else
    sprintf("error:the number of basefunction is incorrect!")

end

end

