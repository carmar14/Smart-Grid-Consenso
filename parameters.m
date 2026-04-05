function P = parameters()

P.l = 6; % DERs

% tiempos (paper sección 4)
P.Ts = 50e-6;
P.tau_ol = 0.5e-3;
P.tau_og = 5e-3;
P.tau_c  = 0.5;
P.Tcomm  = 50e-3;

% horizonte MPC
P.N = 10;

% restricciones
P.u_min = [-0; -1];
P.u_max = [0.8; 1];

% referencias
P.p_ref = 0.15;
P.v_ref = 0.01;

% consenso
P.k0 = 1e-3;

end