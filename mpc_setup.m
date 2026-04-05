function mpc = mpc_setup(Az,Bz,Cz,P)

nx = size(Az,1);
nu = size(Bz,2);

Q = Cz'*diag([5e3*ones(1,nx/2),5e2*ones(1,nx/4)])*Cz;
R = kron(eye(P.l), diag([1e6 1e3]));

mpc.Q = Q;
mpc.R = R;
mpc.N = P.N;

end