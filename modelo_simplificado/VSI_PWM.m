function [m,D,gH,gL,vInv] = VSI_PWM(u,Vdc,t,fsw)
% VSI_PWM
% Genera PWM bipolar para un VSI monofasico de dos niveles.
%
% Entradas:
%   u   = referencia de voltaje del MPC [V]
%   Vdc = tension del bus DC [V]
%   t   = tiempo actual de simulacion [s]
%   fsw = frecuencia de conmutacion [Hz]
%
% Salidas:
%   m    = indice de modulacion [-1,1]
%   D    = duty cycle [0,1]
%   gH   = gate del interruptor superior
%   gL   = gate complementario
%   vInv = tension instantanea ideal del VSI [V]
%
% Modelo PWM bipolar:
%   m = 2*u/Vdc
%   D = (1+m)/2
%
% Portadora triangular:
%   carrier = 1 - 4*abs(frac(t*fsw)-0.5)
%
% Conmutacion:
%   gH = 1 si m >= carrier
%   gL = 1 - gH
%
% Tension de salida:
%   vInv = Vdc*(2*gH-1)
%
% NOTA:
% Este modelo no incluye tiempo muerto. Para una implementacion
% de potencia mas detallada, debe incorporarse dead-time entre
% gH y gL.

% Proteccion contra valores invalidos
if Vdc <= 0 || fsw <= 0
    m = 0;
    D = 0.5;
    gH = 0;
    gL = 0;
    vInv = 0;
    return;
end

% ------------------------------------------------------------
% 1. Indice de modulacion generado por el MPC
% ------------------------------------------------------------
m = 2*u/Vdc;

% Saturacion del indice de modulacion
m = min(max(m,-1),1);

% ------------------------------------------------------------
% 2. Duty cycle equivalente
% ------------------------------------------------------------
D = (1+m)/2;
D = min(max(D,0),1);

% ------------------------------------------------------------
% 3. Portadora triangular normalizada [-1,1]
% ------------------------------------------------------------
phase = t*fsw - floor(t*fsw);
carrier = 1 - 4*abs(phase - 0.5);

% ------------------------------------------------------------
% 4. Generacion de PWM
% ------------------------------------------------------------
if m >= carrier
    gH = 1;
else
    gH = 0;
end

% Interruptor complementario
gL = 1 - gH;

% ------------------------------------------------------------
% 5. Tension instantanea del VSI
% ------------------------------------------------------------
vInv = Vdc*(2*gH - 1);

end
