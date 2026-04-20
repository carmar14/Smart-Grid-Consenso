clear; clc; close all

N=3;
L=5e-3;
R=0.1;
Vg=311;
w0=2*pi*60;

mp=1e-4;
kc=0.002;

Ts=1e-5;
Q=50*eye(2);
Rmpc=0.5*eye(2);

tspan=[0 .2];

Aadj=[0 1 1
      1 0 1
      1 1 0];

alpha=[1;.9;.8];
id_base=[12;8;5];

x0=zeros(N*6,1);

[t,x]=ode45(@(t,x)dyn_dec(...
 t,x,N,L,R,Vg,w0,mp,kc,Aadj,alpha,id_base,Ts,Q,Rmpc),...
 tspan,x0);

%% Potencias nodales
figure; hold on
for i=1:N
 id=x(:,(i-1)*6+1);
 P=Vg*id;
 plot(t,P,'LineWidth',2)
end
grid on
legend('P1','P2','P3')
title('Potencias desacopladas')

%% Corrientes abc
figure
for i=1:N
 id=x(:,(i-1)*6+1);
 iq=x(:,(i-1)*6+2);
 th=x(:,(i-1)*6+5);

 ia=id.*cos(th)-iq.*sin(th);
 ib=id.*cos(th-2*pi/3)-iq.*sin(th-2*pi/3);
 ic=id.*cos(th+2*pi/3)-iq.*sin(th+2*pi/3);

 subplot(N,1,i)
 plot(t,ia,t,ib,t,ic)
 grid on
 title(['Corrientes nodo ',num2str(i)])
end

%% Voltajes por nodo (ya NO PCC)
figure
for i=1:N
 th=x(:,(i-1)*6+5);
 va=Vg*cos(th);
 vb=Vg*cos(th-2*pi/3);
 vc=Vg*cos(th+2*pi/3);

 subplot(N,1,i)
 plot(t,va,t,vb,t,vc)
 grid on
 title(['Voltaje nodo aislado ',num2str(i)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx=dyn_dec(...
~,x,N,L,R,Vg,w0,mp,kc,Aadj,alpha,id_base,Ts,Q,Rmpc)

dx=zeros(size(x));
P=zeros(N,1);

for i=1:N
 idx=(i-1)*6;
 P(i)=Vg*x(idx+1);
end

for i=1:N

 idx=(i-1)*6;
 id=x(idx+1);
 iq=x(idx+2);

 %% consenso SOLO informacional
 s=0;
 for j=1:N
   s=s+Aadj(i,j)*(P(j)/alpha(j)-P(i)/alpha(i));
 end

 id_ref=id_base(i)+kc*s;
 iq_ref=0;

 D=5e-4;

 w=w0-mp*P(i)+kc*s-D*(P(i)-Vg*id_base(i));

 %% MPC LOCAL (independiente por nodo)
 A=[-R/L w;
   -w -R/L];

 B=[1/L 0;
    0 1/L];

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

 Vmax=50;
 vd=max(min(vd,Vmax),-Vmax);
 vq=max(min(vq,Vmax),-Vmax);

 did=(vd-R*id+w*L*iq)/L;
 diq=(vq-R*iq-w*L*id)/L;

 dx(idx+1)=did;
 dx(idx+2)=diq;
 dx(idx+3)=vd;
 dx(idx+4)=vq;
 dx(idx+5)=w;
 dx(idx+6)=0;

end

end