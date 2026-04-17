clear; clc; close all;

%% Parámetros
N = 3;
L = 5e-3;
R = 0.1;
Vg = 311;
w0 = 2*pi*60;

Kp = 20;
Ki = 500;

mp = 0.002;
kc = 0.05;

tau = 0.001; % 10 ms de retardo

A = [0 1 1;
     1 0 1;
     1 1 0];

lags = tau;

tspan = [0 0.3];

%% Historia inicial
history = @(t) zeros(N*6,1);

%% Simulación con retardo
sol = dde23(@(t,x,Z) dynamics_delay(t,x,Z,N,L,R,Vg,w0,Kp,Ki,mp,kc,A), ...
            lags, history, tspan);

t = sol.x';
x = sol.y';

%% Graficar ángulo
figure;
for i=1:N
    theta = x(:,(i-1)*6+5);
    subplot(N,1,i)
    plot(t, theta)
    grid on
    title(['Theta nodo ', num2str(i)])
end

%% Ondas abc
figure;
for i=1:N
    
    id = x(:,(i-1)*6+1);
    iq = x(:,(i-1)*6+2);
    theta = x(:,(i-1)*6+5);
    
    ia = id .* cos(theta) - iq .* sin(theta);
    ib = id .* cos(theta - 2*pi/3) - iq .* sin(theta - 2*pi/3);
    ic = id .* cos(theta + 2*pi/3) - iq .* sin(theta + 2*pi/3);
    
    subplot(N,1,i)
    plot(t, ia, t, ib, t, ic)
    grid on
    title(['abc nodo ', num2str(i)])
end

%% =============================
% Dinámica con retardo
%% =============================
function dx = dynamics_delay(~,x,Z,N,L,R,Vg,w0,Kp,Ki,mp,kc,A)

dx = zeros(size(x));

x_tau = Z(:,1); % estados retardados

P = zeros(N,1);
P_tau = zeros(N,1);

%% Potencia actual
for i=1:N
    idx = (i-1)*6;
    
    id = x(idx+1);
    iq = x(idx+2);
    theta = x(idx+5);
    
    vd = Vg*cos(theta);
    vq = Vg*sin(theta);
    
    P(i) = vd*id + vq*iq;
end

%% Potencia retardada
for i=1:N
    idx = (i-1)*6;
    
    id = x_tau(idx+1);
    iq = x_tau(idx+2);
    theta = x_tau(idx+5);
    
    vd = Vg*cos(theta);
    vq = Vg*sin(theta);
    
    P_tau(i) = vd*id + vq*iq;
end

%% Frecuencia con retardo
omega = zeros(N,1);

for i=1:N
    
    sum_term = 0;
    
    for j=1:N
        sum_term = sum_term + A(i,j)*(P_tau(j) - P(i));
    end
    
    omega(i) = w0 - mp*P(i) + kc*sum_term;
end

%% Dinámica por nodo
for i=1:N
    
    idx = (i-1)*6;
    
    id = x(idx+1);
    iq = x(idx+2);
    xid = x(idx+3);
    xiq = x(idx+4);
    
    w = omega(i);
    
    id_ref = 10;
    iq_ref = 0;
    
    ed = id_ref - id;
    eq = iq_ref - iq;
    
    vd_ctrl = Kp*ed + Ki*xid - w*L*iq;
    vq_ctrl = Kp*eq + Ki*xiq + w*L*id;
    
    did = (vd_ctrl - R*id + w*L*iq)/L;
    diq = (vq_ctrl - R*iq - w*L*id)/L;
    
    dx(idx+1) = did;
    dx(idx+2) = diq;
    dx(idx+3) = ed;
    dx(idx+4) = eq;
    dx(idx+5) = w;
    dx(idx+6) = 0;
end

end