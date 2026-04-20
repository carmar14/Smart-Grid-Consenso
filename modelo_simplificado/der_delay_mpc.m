clc
clear
close all
%% =========================================
% PARAMETROS
%% =========================================

N=3;

L=5e-3;
R=0.1;

Vg=120*sqrt(2);

w0=2*pi*60;

mp=1e-4;

kc=0.002;

Ts=1e-5;

Q=50*eye(2);

Rmpc=0.5*eye(2);

tspan=[0 .2];

%% =========================================
% RED
%% =========================================

Aadj=[0 1 1
      1 0 1
      1 1 0];

alpha=[1
       .9
       .8];

id_base=[12
         8
         5];

%% =========================================
% CONDICION INICIAL
%% =========================================

x0=zeros(N*6,1);

%% =========================================
% SIMULACION
%% =========================================

[t,x]=ode45(...
@(t,x)dynamics_no_delay(...
t,x,...
N,L,R,...
Vg,w0,...
mp,kc,...
Aadj,...
alpha,...
id_base,...
Ts,Q,Rmpc),...
tspan,...
x0);

%% =========================================
% POTENCIAS
%% =========================================

figure
hold on

for i=1:N

id=x(:,(i-1)*6+1);

P=Vg*id;

plot(t,P,'LineWidth',2)

end

grid on

title('Potencias')

legend('P1','P2','P3')


%% =========================================
% POTENCIA TOTAL PCC
%% =========================================

Ptotal=zeros(length(t),1);

for i=1:N

id=x(:,(i-1)*6+1);

Ptotal=Ptotal+Vg*id;

end

figure

plot(t,Ptotal,'LineWidth',2)

grid on

title('Potencia Total PCC')


%% =========================================
% CORRIENTES ABC
%% =========================================

figure

for i=1:N

id=x(:,(i-1)*6+1);

iq=x(:,(i-1)*6+2);

theta=x(:,(i-1)*6+5);

ia=id.*cos(theta)-iq.*sin(theta);

ib=id.*cos(theta-2*pi/3)-iq.*sin(theta-2*pi/3);

ic=id.*cos(theta+2*pi/3)-iq.*sin(theta+2*pi/3);

subplot(N,1,i)

plot(t,ia,t,ib,t,ic)

grid on

title(['Corrientes Nodo ',num2str(i)])

end


%% =========================================
% VOLTAJES RED
%% =========================================

figure

for i=1:N

theta=x(:,(i-1)*6+5);

va=Vg*cos(theta);

vb=Vg*cos(theta-2*pi/3);

vc=Vg*cos(theta+2*pi/3);

subplot(N,1,i)

plot(t,va,t,vb,t,vc)

grid on

title(['Voltaje Red Nodo ',num2str(i)])

end


%% =========================================
% VOLTAJES INVERSORES
%% =========================================

figure

for i=1:N

vd=x(:,(i-1)*6+3);

vq=x(:,(i-1)*6+4);

theta=x(:,(i-1)*6+5);

va=...
vd.*cos(theta)...
-vq.*sin(theta);

vb=...
vd.*cos(theta-2*pi/3)...
-vq.*sin(theta-2*pi/3);

vc=...
vd.*cos(theta+2*pi/3)...
-vq.*sin(theta+2*pi/3);

subplot(N,1,i)

plot(t,va,t,vb,t,vc)

grid on

title(['Voltaje Inversor ',num2str(i)])

end


%% =========================================
% VOLTAJE PCC
%% =========================================

th1=x(:,5);

th2=x(:,11);

th3=x(:,17);


va_pcc=(...
Vg*cos(th1)+...
Vg*cos(th2)+...
Vg*cos(th3))/3;


vb_pcc=(...
Vg*cos(th1-2*pi/3)+...
Vg*cos(th2-2*pi/3)+...
Vg*cos(th3-2*pi/3))/3;


vc_pcc=(...
Vg*cos(th1+2*pi/3)+...
Vg*cos(th2+2*pi/3)+...
Vg*cos(th3+2*pi/3))/3;


figure

plot(t,va_pcc,t,vb_pcc,t,vc_pcc)

grid on

title('Voltaje PCC')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx=dynamics_no_delay(...
~,x,...
N,L,R,...
Vg,w0,...
mp,kc,...
Aadj,...
alpha,...
id_base,...
Ts,Q,Rmpc)

dx=zeros(size(x));

P=zeros(N,1);


%% Potencias

for i=1:N

idx=(i-1)*6;

id=x(idx+1);

P(i)=Vg*id;

end


%% Nodos

for i=1:N

idx=(i-1)*6;

id=x(idx+1);

iq=x(idx+2);


%% consenso

sum_local=0;

for j=1:N

sum_local=...
sum_local...
+Aadj(i,j)*...
(P(j)/alpha(j)...
-P(i)/alpha(i));

end


%% referencias

id_ref=...
id_base(i)...
+kc*sum_local;

iq_ref=0;


%% frecuencia

D=5e-4;

w=...
w0...
-mp*P(i)...
+kc*sum_local...
-D*(P(i)-Vg*id_base(i));


%% MPC

Am=[-R/L w
    -w -R/L];

Bm=[1/L 0
    0 1/L];

Ad=eye(2)+Ts*Am;

Bd=Ts*Bm;


xk=[id;iq];

xref=[id_ref;iq_ref];


H=...
Bd'*Q*Bd...
+Rmpc;

f=...
Bd'*Q*(Ad*xk-xref);


if rcond(H)<1e-10

u=[0;0];

else

u=-H\f;

end


vd=u(1);

vq=u(2);


Vmax=50;

vd=max(min(vd,Vmax),-Vmax);

vq=max(min(vq,Vmax),-Vmax);


%% Planta RL

did=...
(vd...
-R*id...
+w*L*iq)/L;


diq=...
(vq...
-R*iq...
-w*L*id)/L;


dx(idx+1)=did;

dx(idx+2)=diq;

dx(idx+3)=vd;

dx(idx+4)=vq;

dx(idx+5)=w;

dx(idx+6)=0;

end

end