function results = simulate_case_LQR(P,A,B,C,T)

l  = P.l;
Ts = P.Ts;

nx = size(A,1);
ny = size(C,1);
nu = size(B,2);

Nsteps = round(P.Tsim / Ts);

[Ad,Bd] = discretize_model(A,B,Ts);

% LQR
Q = eye(nx);
R = 1e-3*eye(nu);

K = dlqr(Ad,Bd,Q,R);

x = zeros(nx,1);
xe = zeros(ny,1);

y_ref = zeros(ny,1);
for i=1:l
    y_ref(2*i-1) = P.p_ref;
    y_ref(2*i)   = P.v_ref;
end

x_hist = zeros(nx,Nsteps);
y_hist = zeros(ny,Nsteps);
u_hist = zeros(nu,Nsteps);

for k=1:Nsteps

    y = C*x;

    u = -K*x;

    % saturación (como paper)
    u = max(min(u,1),-1);

    x = Ad*x + Bd*u;

    x_hist(:,k) = x;
    y_hist(:,k) = y;
    u_hist(:,k) = u;

end

results.x = x_hist;
results.y = y_hist;
results.u = u_hist;

end