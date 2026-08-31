clc
clear
close all

%----parametros del sistema---
Lf = 1e-3;
Lc = .5e-3;
Cf = 5e-6;

R = 1;
L = 5e-3;

%-----MPC-----
%xp=Ax+Bu+Ed
%y =Cx+Du
%-----|iLf|
%---x=|Vcf|
%-----|ILc|

%---u= Vinv

%--d=vout----perturbacion


A = [0 -1/Lf 0;
     1/Cf 0 -1/Cf;
     0 1/Lc 0];

B = [1/Lf 0 0]';
%---y --_regular corriente hacia la carga y voltaje del capacitor
C = [0 1 0;
     0 0 1];

D = zeros(2,2);

E= [0 0 -1/Lc]';

%%------sistema aumentado para eliminar error en estado estacionario----
% %Aa= [A 0;
%     -Cv 0];
% Ba = [B;0];

%% Sampling Time

Ts = 50e-6;

%% Continuous system

plant_c = ss(A,[B E],C,D);

%% Discrete model

plant = c2d(plant_c,Ts);
%% ============================================================
% MPC OBJECT
% ============================================================

plant.InputGroup.MV = 1;
plant.InputGroup.MD = 2;

PredictionHorizon = 20; %Np

ControlHorizon = 5;   %Nu

mpcobj = mpc(plant,Ts,...
             PredictionHorizon,...
             ControlHorizon);

%% ============================================================
% WEIGHTS
% ============================================================

mpcobj.Weights.OutputVariables = [10 1];

mpcobj.Weights.ManipulatedVariables = 0;

mpcobj.Weights.ManipulatedVariablesRate = 0.1;

%% ============================================================
% INPUT CONSTRAINTS
% ============================================================
Vdc = 400;
mpcobj.MV.Min = -Vdc/2;

mpcobj.MV.Max = Vdc/2;

mpcobj.MV.RateMin = -100;

mpcobj.MV.RateMax = 100;

%% ============================================================
% OUTPUT CONSTRAINTS (optional)
% ============================================================
Vp = sqrt(2)*120;
% voltage limitation
mpcobj.OV(1).Min = -Vp;
mpcobj.OV(1).Max = Vp;

% current limitation
Imax = 20;
mpcobj.OV(2).Min = -Imax;
mpcobj.OV(2).Max = Imax;

%% ============================================================
% NOMINAL CONDITIONS
% ============================================================

% mpcobj.Model.Nominal.U = 0;
% 
% mpcobj.Model.Nominal.Y = 0;
% 
% mpcobj.Model.Nominal.X = zeros(size(plant.A,1),1);

%% ============================================================
% ESTIMATOR
% ============================================================

setEstimator(mpcobj,'default');
setoutdist(mpcobj,'integrators');

%% ============================================================
% SAVE OBJECT
% ============================================================

save MPC_Local_Controller mpcobj

disp('MPC created successfully')