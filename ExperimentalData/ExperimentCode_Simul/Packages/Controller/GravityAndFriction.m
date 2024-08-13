function uff = GravityAndFriction(q,dq)
jiaodu = q.';
sudu = dq.';

guanjieliju_1 = -51.6115*sin(jiaodu(1)) -21.2854*sin(jiaodu(1)+jiaodu(2)) ...
                -6.0632*sin(jiaodu(1)+jiaodu(2)+jiaodu(3))...
                +27.7593*satr(sudu(1))*0.0+ 14.7012*sudu(1)*1.0;
          
guanjieliju_2 = -23.1921*sin(jiaodu(1)+jiaodu(2)) - 5.4261*sin(jiaodu(1)+jiaodu(2)+jiaodu(3))...
                +7.8431*satr(sudu(2))*0.0 + 11.5956*sudu(2)*1.0;%  

guanjieliju_3 = 5.7248*sin(jiaodu(1)+jiaodu(2)+jiaodu(3))...
                -2.5745*satr(sudu(3))*0.0 - 3.7623*sudu(3)*1.0;

uff = [guanjieliju_1; guanjieliju_2; guanjieliju_3];
end

