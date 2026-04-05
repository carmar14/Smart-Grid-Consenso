function u = mpc_solve(z,mpc,umin,umax)

H = blkdiag(kron(eye(mpc.N),mpc.R));
f = zeros(size(H,1),1);

Aineq = [];
bineq = [];

lb = repmat(umin,mpc.N,1);
ub = repmat(umax,mpc.N,1);

U = quadprog(H,f,Aineq,bineq,[],[],lb,ub);

u = U(1:length(umin));

end