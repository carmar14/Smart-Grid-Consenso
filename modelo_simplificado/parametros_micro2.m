clc
clear
close all

%% ============================================================
% PARAMETROS DEL SISTEMA
% ============================================================

Lf = 1e-3;
Lc = 0.5e-3;
Cf = 5e-6;

% Carga RL
R = 1;
L = 5e-3;

% Tiempo de muestreo
Ts = 50e-6;

% Bus DC
Vdc = 400;

%% ============================================================
% PARAMETROS DEL MPC
% ============================================================

PredictionHorizon = 20;
ControlHorizon = 5;

%% ============================================================
% MODELO LOCAL DE CADA INVERSOR
% ============================================================

A = [0 -1/Lf 0;
     1/Cf 0 -1/Cf;
     0 1/Lc 0];

B = [1/Lf;
     0;
     0];

E = [0;
     0;
     -1/Lc];

%% ============================================================
% ECUACION DEL PCC
%
% vPCC =
%
% R*(Ic1+Ic2+Ic3)
% +
% (L/Lc)*(Vc1+Vc2+Vc3)
%
% ---------------------------------
%
% 1 + 3*(L/Lc)
% ============================================================

alpha = L/Lc;

Kpcc = 1/(1 + 3*alpha);

%% ============================================================
% PARTE LOCAL DEL PCC
% ============================================================

C_local = Kpcc*[0 alpha R];

%% ============================================================
% MODELO LOCAL REFORMULADO
%
% dx/dt =
%
% (A + E*C_local)x
% + B*u
% + E*d
%
% vPCC = C_local*x + d
% ============================================================

A_local = A + E*C_local;

B_local = B;

E_local = E;

D_local = [0 1];

%% ============================================================
% CREAR MPC 1
% ============================================================

plant_c1 = ss(...
    A_local,...
    [B_local E_local],...
    C_local,...
    D_local);

plant1 = c2d(plant_c1,Ts);

plant1 = setmpcsignals(...
    plant1,...
    'MV',1,...
    'MD',2);

mpc1 = mpc(...
    plant1,...
    Ts,...
    PredictionHorizon,...
    ControlHorizon);

%% ============================================================
% CREAR MPC 2
% ============================================================

plant_c2 = ss(...
    A_local,...
    [B_local E_local],...
    C_local,...
    D_local);

plant2 = c2d(plant_c2,Ts);

plant2 = setmpcsignals(...
    plant2,...
    'MV',1,...
    'MD',2);

mpc2 = mpc(...
    plant2,...
    Ts,...
    PredictionHorizon,...
    ControlHorizon);

%% ============================================================
% CREAR MPC 3
% ============================================================

plant_c3 = ss(...
    A_local,...
    [B_local E_local],...
    C_local,...
    D_local);

plant3 = c2d(plant_c3,Ts);

plant3 = setmpcsignals(...
    plant3,...
    'MV',1,...
    'MD',2);

mpc3 = mpc(...
    plant3,...
    Ts,...
    PredictionHorizon,...
    ControlHorizon);

%% ============================================================
% CONFIGURACION DE LOS TRES MPC
% ============================================================

% ------------------------------------------------------------
% PESOS
% ------------------------------------------------------------

mpc1.Weights.OutputVariables = 10;
mpc2.Weights.OutputVariables = 10;
mpc3.Weights.OutputVariables = 10;

mpc1.Weights.ManipulatedVariables = 1;
mpc2.Weights.ManipulatedVariables = 1;
mpc3.Weights.ManipulatedVariables = 1;

mpc1.Weights.ManipulatedVariablesRate = 0.1;
mpc2.Weights.ManipulatedVariablesRate = 0.1;
mpc3.Weights.ManipulatedVariablesRate = 0.1;

% ------------------------------------------------------------
% RESTRICCIONES MV
% ------------------------------------------------------------
mpc1.MV.Min = -Vdc/2;
mpc1.MV.Max =  Vdc/2;
mpc1.MV.RateMin = -100;
mpc1.MV.RateMax =  100;

mpc2.MV.Min = -Vdc/2;
mpc2.MV.Max =  Vdc/2;
mpc2.MV.RateMin = -100;
mpc2.MV.RateMax =  100;

mpc3.MV.Min = -Vdc/2;
mpc3.MV.Max =  Vdc/2;
mpc3.MV.RateMin = -100;
mpc3.MV.RateMax =  100;


% ------------------------------------------------------------
% RESTRICCIONES PCC
% ------------------------------------------------------------

Vpcc_peak = sqrt(2)*120;

mpc1.OV.Min = -Vpcc_peak;
mpc1.OV.Max =  Vpcc_peak;

mpc2.OV.Min = -Vpcc_peak;
mpc2.OV.Max =  Vpcc_peak;

mpc3.OV.Min = -Vpcc_peak;
mpc3.OV.Max =  Vpcc_peak;

%% ============================================================
% ESTIMADORES
% ============================================================

setEstimator(mpc1,'default');
setEstimator(mpc2,'default');
setEstimator(mpc3,'default');

%% ============================================================
% MODELOS DE PERTURBACION
% ============================================================

setoutdist(mpc1,'integrators');
setoutdist(mpc2,'integrators');
setoutdist(mpc3,'integrators');

%% ============================================================
% REVISAR MPC
% ============================================================

disp('============================================================')
disp('REVISION MPC1')
disp('============================================================')

%review(mpc1)

disp('============================================================')
disp('REVISION MPC2')
disp('============================================================')

%review(mpc2)

disp('============================================================')
disp('REVISION MPC3')
disp('============================================================')

%review(mpc3)

%% ============================================================
% SIMULACION
% ============================================================

Tsim = 0.2;

N = round(Tsim/Ts);

t = (0:N-1)'*Ts;

%% ============================================================
% REFERENCIA PCC
% ============================================================

f = 60;

Vpcc_peak = sqrt(2)*120;

r = Vpcc_peak*sin(2*pi*f*t);

%% ============================================================
% ESTADOS INICIALES
%
% x1 = [if1 Vc1 Ic1]'
% x2 = [if2 Vc2 Ic2]'
% x3 = [if3 Vc3]'
% ============================================================

x1 = zeros(3,1);
x2 = zeros(3,1);
x3 = zeros(3,1);

%% ============================================================
% ESTADOS INTERNOS DE LOS MPC
% ============================================================

state1 = mpcstate(mpc1);
state2 = mpcstate(mpc2);
state3 = mpcstate(mpc3);

%% ============================================================
% ALMACENAMIENTO
% ============================================================

x1_hist = zeros(3,N);
x2_hist = zeros(3,N);
x3_hist = zeros(3,N);

vPCC_hist = zeros(1,N);

d1_hist = zeros(1,N);
d2_hist = zeros(1,N);
d3_hist = zeros(1,N);

u1_hist = zeros(1,N);
u2_hist = zeros(1,N);
u3_hist = zeros(1,N);

%% ============================================================
% LAZO CERRADO
% ============================================================

for k = 1:N

    %% ========================================================
    % 1. CALCULAR VOLTAJE PCC REAL
    % ========================================================

    Vc1 = x1(2);
    Ic1 = x1(3);

    Vc2 = x2(2);
    Ic2 = x2(3);

    Vc3 = x3(2);
    Ic3 = x3(3);

    vPCC = Kpcc*(...
        alpha*(Vc1 + Vc2 + Vc3) + ...
        R*(Ic1 + Ic2 + Ic3));

    %% ========================================================
    % 2. CALCULAR PERTURBACIONES LOCALES
    % ========================================================

    d1 = Kpcc*(...
        alpha*(Vc2 + Vc3) + ...
        R*(Ic2 + Ic3));

    d2 = Kpcc*(...
        alpha*(Vc1 + Vc3) + ...
        R*(Ic1 + Ic3));

    d3 = Kpcc*(...
        alpha*(Vc1 + Vc2) + ...
        R*(Ic1 + Ic2));

    %% ========================================================
    % 3. MEDICION PARA CADA MPC
    % ========================================================

    y1 = vPCC;
    y2 = vPCC;
    y3 = vPCC;

    %% ========================================================
    % 4. MPC1
    % ========================================================

    u1 = mpcmove(...
        mpc1,...
        state1,...
        y1,...
        r(k),...
        d1);

    %% ========================================================
    % 5. MPC2
    % ========================================================

    u2 = mpcmove(...
        mpc2,...
        state2,...
        y2,...
        r(k),...
        d2);

    %% ========================================================
    % 6. MPC3
    % ========================================================

    u3 = mpcmove(...
        mpc3,...
        state3,...
        y3,...
        r(k),...
        d3);

    %% ========================================================
    % 7. SATURACION
    % ========================================================

    u1 = min(max(u1,-Vdc/2),Vdc/2);
    u2 = min(max(u2,-Vdc/2),Vdc/2);
    u3 = min(max(u3,-Vdc/2),Vdc/2);

    %% ========================================================
    % 8. PLANTA FISICA
    %
    % IMPORTANTE:
    %
    % Aqui usamos el vPCC REAL.
    %
    % dx/dt = A*x + B*u + E*vPCC
    % ========================================================

    x1_next = plant1.A*x1 + ...
              plant1.B(:,1)*u1 + ...
              plant1.B(:,2)*vPCC;

    x2_next = plant2.A*x2 + ...
              plant2.B(:,1)*u2 + ...
              plant2.B(:,2)*vPCC;

    x3_next = plant3.A*x3 + ...
              plant3.B(:,1)*u3 + ...
              plant3.B(:,2)*vPCC;

    %% ========================================================
    % 9. GUARDAR DATOS
    % ========================================================

    x1_hist(:,k) = x1;
    x2_hist(:,k) = x2;
    x3_hist(:,k) = x3;

    vPCC_hist(k) = vPCC;

    d1_hist(k) = d1;
    d2_hist(k) = d2;
    d3_hist(k) = d3;

    u1_hist(k) = u1;
    u2_hist(k) = u2;
    u3_hist(k) = u3;

    %% ========================================================
    % 10. ACTUALIZAR ESTADOS
    % ========================================================

    x1 = x1_next;
    x2 = x2_next;
    x3 = x3_next;

end

%% ============================================================
% ERROR DE SEGUIMIENTO
% ============================================================

error = r' - vPCC_hist;

%% ============================================================
% RMS DEL PCC
% ============================================================

idx = round(N/2):N;

Vrms_PCC = rms(vPCC_hist(idx));

disp('============================================================')
disp('RESULTADOS')
disp('============================================================')

fprintf('Voltaje PCC RMS = %.4f V\n',Vrms_PCC);

fprintf('Referencia RMS = %.4f V\n',120);

fprintf('Error RMS = %.4f V\n',Vrms_PCC - 120);

%% ============================================================
% GRAFICA 1
% VOLTAJE PCC
% ============================================================

figure

plot(t,vPCC_hist,'LineWidth',1.5)

hold on

plot(t,r,'--','LineWidth',1.2)

grid on

xlabel('Time [s]')
ylabel('v_{PCC} [V]')

legend('v_{PCC}','v_{PCC}^{*}')

title('Lazo Cerrado - Voltaje PCC')

%% ============================================================
% GRAFICA 2
% ERROR DE SEGUIMIENTO
% ============================================================

figure

plot(t,error,'LineWidth',1.5)

grid on

xlabel('Time [s]')
ylabel('Error [V]')

title('Error de Seguimiento del PCC')

%% ============================================================
% GRAFICA 3
% TENSIONES DE LOS INVERSORES
% ============================================================

figure

plot(t,x1_hist(2,:),'LineWidth',1.5)

hold on

plot(t,x2_hist(2,:),'LineWidth',1.5)

plot(t,x3_hist(2,:),'LineWidth',1.5)

grid on

xlabel('Time [s]')
ylabel('V_C [V]')

legend('V_{C1}','V_{C2}','V_{C3}')

title('Tensiones de los Capacitores')

%% ============================================================
% GRAFICA 4
% CORRIENTES DE LOS INVERSORES
% ============================================================

figure

plot(t,x1_hist(3,:),'LineWidth',1.5)

hold on

plot(t,x2_hist(3,:),'LineWidth',1.5)

plot(t,x3_hist(3,:),'LineWidth',1.5)

grid on

xlabel('Time [s]')
ylabel('I_c [A]')

legend('I_{c1}','I_{c2}','I_{c3}')

title('Corrientes de los Tres Inversores')

%% ============================================================
% GRAFICA 5
% ACCIONES DE CONTROL
% ============================================================

figure

plot(t,u1_hist,'LineWidth',1.5)

hold on

plot(t,u2_hist,'LineWidth',1.5)

plot(t,u3_hist,'LineWidth',1.5)

grid on

xlabel('Time [s]')
ylabel('V_{inv} [V]')

legend('V_{inv1}','V_{inv2}','V_{inv3}')

title('Acciones de Control de los Tres MPC')

%% ============================================================
% GRAFICA 6
% PERTURBACIONES LOCALES
% ============================================================

figure

plot(t,d1_hist,'LineWidth',1.5)

hold on

plot(t,d2_hist,'LineWidth',1.5)

plot(t,d3_hist,'LineWidth',1.5)

grid on

xlabel('Time [s]')
ylabel('d_i [V]')

legend('d_1','d_2','d_3')

title('Perturbaciones Locales')