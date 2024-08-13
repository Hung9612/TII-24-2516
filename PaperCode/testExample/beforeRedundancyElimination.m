clear
clc
%%
syms q1 q2 q3;
syms dq1 dq2 dq3;
syms ddq1 ddq2 ddq3 L1 L2 L3 L4 g;
n = 3;
q = [q1,q2,q3].';
dq = [dq1 dq2 dq3].';
ddq = [ddq1 ddq2 ddq3].';
L = [L1 L2 L3 L4].';
%%
e1 = [xrot(q1),zeros(1,3)';0,0,0,1];
e2 = [zrot(q2),zeros(1,3)';0,0,0,1];
e3 = [xrot(q3),zeros(1,3)';0,0,0,1];

T0 = transl([0,0,0.3245]);
T1 = transl([0,0,0.35]);
T2 = transl([0,0,0.35]);
T3 = transl([0,0,0.275]);
T01 = T0*e1;
T02 = T01*T1*e2;
T03 = T02*T2*T3*e3;
T = {T01,T02,T03};
%% 
D_mtrx = zeros(3,3)*q1;
for i = 1:n
    for j = 1:n
        p = max(i,j);
        d_ij = q1*0;
        
        for k = p:n
            
            partial_Tpi = diff(T{k}, q(i));
            partial_Tpj = diff(T{k}.', q(j));

            J_p = fetch_para_before(k, 3);

            out = trace(partial_Tpi*J_p*partial_Tpj);
            
            d_ij = d_ij + out;%Tr_dij(k,q(i),q(j),T{k});
        end
        
        D_mtrx(i,j) = d_ij;
    end
end
tauD = simplify(D_mtrx*ddq,'step',1);

TauD=expand(tauD);
%% 


C_mtrx = zeros(3,3)*q1;
for i = 1:n
    C_jk = zeros(3,3)*q1;
    for j = 1:n
        for k = 1:n
            c_ijk = q1*0;
            p = max([i,j,k]);
             for t = p:n
                 
                  partial_Tpi = diff(T{t}.', q(i));
    
                  partial_Tpj = diff(T{t}, q(j));
                  partial_Tpjk = diff(partial_Tpj,q(k));

                  J_p = fetch_para(t, 3);

                  out = trace(partial_Tpjk*J_p*partial_Tpi);
                  c_ijk = c_ijk + out;% Tr_cijk(t,q(i),q(j),q(k),T{t});
             end
             
             C_jk(j,k) = c_ijk;
        end
    end
    C_mtrx(i,:) = dq.'*C_jk;
end
tau_C = C_mtrx*dq;
tau_C = simplify(tau_C,'step',1);
TauC = expand(tau_C);
%% 
g = [0,0,-9.81,0].';

G_vec = q*0;
for i = 1:n
    temp = q1*0;
    for p = i:n
        [m_p,rc_p] = cal_Gpara_before(p,3);
        partial_p = diff(T{p},q(i));
        
        temp = temp + g.' * partial_p * rc_p;
    end
    G_vec(i) = temp;
end
tau_G = simplify(G_vec,'step',1);
TauG = expand(tau_G);

%%

tau = D_mtrx*ddq+ C_mtrx*dq + G_vec;
Tau = simplify(tau,'step',1);
