function [m,D] = VSI_PWM_Modulator(u,Vdc)
% VSI_PWM_Modulator
% Convierte la referencia de voltaje del MPC en indice de
% modulacion y duty cycle para la etapa VSI/PWM.
%
% Entradas:
%   u   = referencia de voltaje del inversor [V]
%   Vdc = voltaje del bus DC [V]
%
% Salidas:
%   m = indice de modulacion [-1,1]
%   D = duty cycle [0,1]
%
% Relaciones:
%   m = 2*u/Vdc
%   D = (1+m)/2
%
% Para el modelo promedio:
%   v_inv_avg = (Vdc/2)*m = u

% Proteccion contra Vdc no valido
if Vdc <= 0
    m = 0;
    D = 0.5;
    return;
end

% Indice de modulacion
m = 2*u/Vdc;

% Saturacion del indice de modulacion
m = min(max(m,-1),1);

% Duty cycle equivalente
D = (1+m)/2;

% Saturacion del duty cycle
D = min(max(D,0),1);

end
