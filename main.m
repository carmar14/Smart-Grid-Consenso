clear; clc; close all;

addpath(genpath(pwd));

% =========================
% PARÁMETROS
% =========================
P = parameters();
P.Tsim = 2; % tiempo total simulación

l = P.l;

% =========================
% MODELOS
% =========================
[Aj,Bj,Cj] = build_DER_exact();
[A,B,Ctilde] = build_global_model(Aj,Bj,Cj,l);

% matriz de acoplamiento (paper)
T = build_T_exact();

% salida completa
C = T * Ctilde;

% =========================
% REFERENCIA
% =========================
y_ref = zeros(2*l,1);
for i=1:l
    y_ref(2*i-1) = P.p_ref;
    y_ref(2*i)   = P.v_ref;
end

% =========================
% CASE I
% Nominal sin comunicación
% =========================
disp('Running Case I...');

P_case = P;
P_case.Tcomm = Inf; % sin sample & hold

res1 = simulate_case(P_case,A,B,C,T);

figure;
plot(res1.y');
title('Case I - Outputs');
xlabel('Time step'); ylabel('y');
legend_strings(l);

% =========================
% CASE II
% Con comunicación
% =========================
disp('Running Case II...');

P_case = P; % Tcomm activo

res2 = simulate_case(P_case,A,B,C,T);

figure;
plot(res2.y');
title('Case II - With Communication');
xlabel('Time step'); ylabel('y');
legend_strings(l);

% =========================
% CASE III
% MPC -> LQR
% =========================
disp('Running Case III...');

P_case = P;
P_case.use_LQR = true;

res3 = simulate_case_LQR(P_case,A,B,C,T);

figure;
plot(res3.y');
title('Case III - LQR Control');
xlabel('Time step'); ylabel('y');
legend_strings(l);

% =========================
% CASE IV
% Parámetros perturbados
% =========================
disp('Running Case IV...');

T_pert = build_T_perturbed(); % matriz T' del paper
C_pert = T_pert * Ctilde;

res4 = simulate_case(P,A,B,C_pert,T_pert);

figure;
plot(res4.y');
title('Case IV - Off-nominal');
xlabel('Time step'); ylabel('y');
legend_strings(l);

% =========================
% COMPARACIÓN GLOBAL
% =========================
figure;

subplot(2,2,1)
plot(res1.y'); title('Case I');

subplot(2,2,2)
plot(res2.y'); title('Case II');

subplot(2,2,3)
plot(res3.y'); title('Case III');

subplot(2,2,4)
plot(res4.y'); title('Case IV');

sgtitle('Comparison of All Cases');

disp('Simulation complete.');