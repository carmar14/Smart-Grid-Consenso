clear; clc; close all;

%% =============================
% Parámetros
%% =============================
N = 3;                  % número de nodos
L = 5e-3;
R = 0.1;
Vg = 120*sqrt(2);               % voltaje pico (120 RMS)
w0 = 2*pi*60;

Kp = 20; %----Kp mas pequeño genera sobreimpulsos
Ki = 500;

% Droop
mp = 0.002;

% Consenso
kc = 0.05;

% Matriz de adyacencia (todos conectados)
A = [0 1 1;
     1 0 1;
     1 1 0];

%% =============================
% Condiciones iniciales
%% =============================
x0 = zeros(N*6,1);

% Pequeña perturbación inicial en ángulo
for i=1:N
    x0((i-1)*6 + 5) = rand*0.1;
end

tspan = [0 0.3];

%% =============================
% Simulación
%% =============================
[t,x] = ode45(@(t,x) dynamics(t,x,N,L,R,Vg,w0,Kp,Ki,mp,kc,A), tspan, x0);

%% =============================
% Graficar ángulos
%% =============================
figure;
for i=1:N
    theta = x(:,(i-1)*6+5);
    subplot(N,1,i)
    plot(t, theta, 'LineWidth',1.5)
    grid on
    title(['Ángulo θ nodo ', num2str(i)])
end

%% =============================
% Graficar corrientes dq
%% =============================
figure;
for i=1:N
    id = x(:,(i-1)*6+1);
    iq = x(:,(i-1)*6+2);
    
    subplot(N,1,i)
    plot(t, id, 'LineWidth',1.5); hold on;
    plot(t, iq, '--', 'LineWidth',1.5);
    grid on
    legend('i_d','i_q')
    title(['Corrientes dq nodo ', num2str(i)])
end

%% =============================
% Reconstrucción abc (ONDAS SENO)
%% =============================
figure;

for i = 1:N
    
    id = x(:,(i-1)*6+1);
    iq = x(:,(i-1)*6+2);
    theta = x(:,(i-1)*6+5);
    
    ia = id .* cos(theta) - iq .* sin(theta);
    ib = id .* cos(theta - 2*pi/3) - iq .* sin(theta - 2*pi/3);
    ic = id .* cos(theta + 2*pi/3) - iq .* sin(theta + 2*pi/3);
    
    subplot(N,1,i)
    plot(t, ia, 'LineWidth',1.2); hold on
    plot(t, ib, 'LineWidth',1.2);
    plot(t, ic, 'LineWidth',1.2);
    grid on
    title(['Corrientes trifásicas nodo ', num2str(i)])
end

figure;

for i = 1:N
    
    theta = x(:,(i-1)*6+5);
    
    Va = Vg * cos(theta);
    Vb = Vg * cos(theta - 2*pi/3);
    Vc = Vg * cos(theta + 2*pi/3);
    
    subplot(N,1,i)
    plot(t, Va, t, Vb, t, Vc)
    grid on
    title(['Voltajes nodo ', num2str(i)])
end


%% =============================
% Dinámica del sistema
%% =============================
function dx = dynamics(~,x,N,L,R,Vg,w0,Kp,Ki,mp,kc,A)

dx = zeros(size(x));

omega = zeros(N,1);
P = zeros(N,1);

%% =============================
% Calcular potencia por nodo
%% =============================
for i=1:N
    
    idx = (i-1)*6;
    
    id = x(idx+1);
    iq = x(idx+2);
    theta = x(idx+5);
    
    % Voltaje ficticio de red
    vd = Vg*cos(theta);
    vq = Vg*sin(theta);
    
    P(i) = vd*id + vq*iq;
end

%% =============================
% Frecuencia con droop + consenso
%% =============================
for i=1:N
    
    sum_term = 0;
    
    for j=1:N
        sum_term = sum_term + A(i,j)*(P(j) - P(i));
    end
    
    omega(i) = w0 - mp*P(i) + kc*sum_term;
end

%% =============================
% Dinámica por nodo
%% =============================
for i=1:N
    
    idx = (i-1)*6;
    
    id = x(idx+1);
    iq = x(idx+2);
    xid = x(idx+3);
    xiq = x(idx+4);
    theta = x(idx+5);
    
    w = omega(i);
    
    % Referencias
    id_ref = 7;
    iq_ref = 0;
    
    % Errores
    ed = id_ref - id;
    eq = iq_ref - iq;
    
    % Control PI + desacoplo
    vd_ctrl = Kp*ed + Ki*xid - w*L*iq;
    vq_ctrl = Kp*eq + Ki*xiq + w*L*id;
    
    % Dinámica RL en dq
    did = (vd_ctrl - R*id + w*L*iq)/L;
    diq = (vq_ctrl - R*iq - w*L*id)/L;
    
    % Estados
    dx(idx+1) = did;
    dx(idx+2) = diq;
    dx(idx+3) = ed;
    dx(idx+4) = eq;
    
    % Ángulo
    dx(idx+5) = w;
    
    % Estado auxiliar
    dx(idx+6) = 0;
end

end