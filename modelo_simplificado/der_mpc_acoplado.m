clear;clc;close all

N=3;

%% Inversores
L=5e-3;
R=0.1;
w0=2*pi*60;

Vnom=311;

mp=1e-4;
kc=0.002;

Ts=1e-5;
Q=50*eye(2);
Rmpc=.5*eye(2);

%% Lineas
Rl=0.05;
Ll=2e-3;

%% Red completa
Aadj=[0 1 1
      1 0 1
      1 1 0];

alpha=[1;.9;.8];
id_base=[12;8;5];

%% 18 estados inversores + 6 lineas
x0=zeros(24,1);

tspan=[0 .2];

[t,x]=ode45(@(t,x)dyn_grid(...
 t,x,N,L,R,Rl,Ll,Vnom,w0,...
 mp,kc,Aadj,alpha,id_base,Ts,Q,Rmpc),...
 tspan,x0);

%% Potencias
figure; hold on
for i=1:N
 id=x(:,(i-1)*6+1);
 iq=x(:,(i-1)*6+2);
 vd=x(:,(i-1)*6+3);
 vq=x(:,(i-1)*6+4);

 P=vd.*id+vq.*iq;
 plot(t,P,'LineWidth',2)
end
grid on
legend('P1','P2','P3')
title('Potencias con red acoplada')

%% Corrientes lineas
figure
plot(t,x(:,19),t,x(:,21),t,x(:,23))
grid on
title('Corrientes d en lineas')
legend('i12','i23','i31')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx=dyn_grid(...
~,x,N,L,R,Rl,Ll,Vnom,w0,...
mp,kc,Aadj,alpha,id_base,Ts,Q,Rmpc)

dx=zeros(24,1);
P=zeros(N,1);

%% estados lineas
id12=x(19); iq12=x(20);
id23=x(21); iq23=x(22);
id31=x(23); iq31=x(24);

%% potencias locales
for i=1:N
 idx=(i-1)*6;
 id=x(idx+1);
 iq=x(idx+2);
 vd=x(idx+3);
 vq=x(idx+4);

 P(i)=vd*id+vq*iq;
end

%% inversores
for i=1:N

idx=(i-1)*6;
id=x(idx+1);
iq=x(idx+2);

s=0;
for j=1:N
s=s+Aadj(i,j)*(P(j)/alpha(j)-P(i)/alpha(i));
end

id_ref=id_base(i)+kc*s;
iq_ref=0;

w=w0-mp*P(i)+kc*s;

A=[-R/L w;-w -R/L];
B=[1/L 0;0 1/L];

Ad=eye(2)+Ts*A;
Bd=Ts*B;

xk=[id;iq];
xr=[id_ref;iq_ref];

H=Bd'*Q*Bd+Rmpc;
f=Bd'*Q*(Ad*xk-xr);

if rcond(H)<1e-10
u=[0;0];
else
u=-H\f;
end

vd=u(1);
vq=u(2);

%% acoplamiento por lineas
if i==1
 iinj_d=id12-id31;
 iinj_q=iq12-iq31;
elseif i==2
 iinj_d=id23-id12;
 iinj_q=iq23-iq12;
else
 iinj_d=id31-id23;
 iinj_q=iq31-iq23;
end

did=(vd-R*id+w*L*iq-iinj_d)/L;
diq=(vq-R*iq-w*L*id-iinj_q)/L;

dx(idx+1)=did;
dx(idx+2)=diq;
dx(idx+3)=vd;
dx(idx+4)=vq;
dx(idx+5)=w;
dx(idx+6)=0;

end

%% voltajes nodales (aprox desde estados control)
v1d=x(3);  v1q=x(4);
v2d=x(9);  v2q=x(10);
v3d=x(15); v3q=x(16);

w=w0;

%% linea 12

dx(19)=(v1d-v2d-Rl*id12+w*Ll*iq12)/Ll;
dx(20)=(v1q-v2q-Rl*iq12-w*Ll*id12)/Ll;

%% linea 23

dx(21)=(v2d-v3d-Rl*id23+w*Ll*iq23)/Ll;
dx(22)=(v2q-v3q-Rl*iq23-w*Ll*id23)/Ll;

%% linea 31

dx(23)=(v3d-v1d-Rl*id31+w*Ll*iq31)/Ll;
dx(24)=(v3q-v1q-Rl*iq31-w*Ll*id31)/Ll;

end
