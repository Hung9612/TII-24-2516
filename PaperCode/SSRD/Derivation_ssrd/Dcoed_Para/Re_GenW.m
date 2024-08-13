function  Re_GenW(W,Q,dQ,ddQ)
Wmatrix =W;

matlabFunction(Wmatrix,'File','W_tau','Var',{Q,dQ,ddQ});
end

