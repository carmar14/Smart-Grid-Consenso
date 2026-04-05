function results = simulate_case(P,A,B,C,T)

l  = P.l;
Ts = P.Ts;

nx = size(A,1);     % 6l
ny = size(C,1);     % 2l
nu = size(B,2);     % 2l

Nsteps = round(P.Tsim / Ts);

% Discretización
[Ad,Bd] = discretize_model(A,B,Ts);

% Modelo extendido
[Az,Bz,Cz] = build_extended_model(Ad,Bd,C);

% MPC
mpc = mpc_setup(Az,Bz,Cz,P);

% =========================
% INICIALIZACIÓN
% =========================
x  = zeros(nx,1);
xe = zeros(ny,1);

% Observadores locales (uno por DER)
xhat_local = zeros(nx,l);

% Observadores consenso (cada DER estima todo el estado)
xhat_global = zeros(nx,l);

% comunicación (sample & hold)
xhat_comm = xhat_global;

comm_counter = 0;
comm_steps = round(P.Tcomm / Ts);

% referencia
y_ref = zeros(ny,1);
for i=1:l
    y_ref(2*i-1) = P.p_ref;
    y_ref(2*i)   = P.v_ref;
end

% almacenamiento
x_hist = zeros(nx,Nsteps);
y_hist = zeros(ny,Nsteps);
u_hist = zeros(nu,Nsteps);

% =========================
% MATRIZ DE VECINOS (FIG 3)
% =========================
neighbors = {
    [2 3 4],
    [1 3],
    [1 2],
    [1 5 6],
    [4 6],
    [4 5]
};

% =========================
% LOOP DE SIMULACIÓN
% =========================
for k = 1:Nsteps

    % salida real
    y = C*x;

    % =====================
    % OBSERVADORES LOCALES
    % =====================
    for i = 1:l
        idx = (i-1)*6+1:i*6;
        y_idx = (i-1)*2+1:i*2;

        xi_hat = xhat_local(idx,i);
        ui = zeros(2,1); % se actualizará después

        Ci = C(y_idx,idx);

        % Ganancia observador (rápida)
        Ki = place(Ad(idx,idx)',Ci', ...
             exp(-Ts/P.tau_ol)*ones(6,1))';

        xhat_local(idx,i) = ...
            Ad(idx,idx)*xi_hat + Bd(idx,:)*zeros(nu,1) + ...
            Ki*(y(y_idx) - Ci*xi_hat);
    end

    % =====================
    % CONSENSO
    % =====================
    for i = 1:l

        xi = xhat_global(:,i);

        sum_term = zeros(nx,1);

        for j = neighbors{i}
            sum_term = sum_term + ...
                (xhat_comm(:,j) - xhat_comm(:,i));
        end

        % matriz selección Π_i
        Pi = zeros(nx);
        idx = (i-1)*6+1:i*6;
        Pi(idx,idx) = eye(6);

        local_term = Pi*(xhat_local(:,i) - xi);

        xhat_global(:,i) = ...
            Ad*xi + Bd*zeros(nu,1) + ...
            P.k0*(sum_term + local_term);
    end

    % =====================
    % SAMPLE & HOLD
    % =====================
    if comm_counter == 0
        xhat_comm = xhat_global;
    end

    comm_counter = comm_counter + 1;
    if comm_counter >= comm_steps
        comm_counter = 0;
    end

    % =====================
    % CONTROL MPC (cada DER)
    % =====================
    u = zeros(nu,1);

    for i = 1:l

        z = [xhat_global(:,i); xe];

        u_i = mpc_solve(z,mpc,P.u_min,P.u_max);

        % cada DER aplica SOLO su parte
        idx_u = (i-1)*2+1:i*2;
        u(idx_u) = u_i(idx_u);
    end

    % =====================
    % DINÁMICA REAL
    % =====================
    x = Ad*x + Bd*u;

    % =====================
    % INTEGRADOR (MPC)
    % =====================
    xe = xe + (y_ref - y)*Ts;

    % =====================
    % GUARDAR
    % =====================
    x_hist(:,k) = x;
    y_hist(:,k) = y;
    u_hist(:,k) = u;

end

% =========================
% OUTPUT
% =========================
results.x = x_hist;
results.y = y_hist;
results.u = u_hist;

end