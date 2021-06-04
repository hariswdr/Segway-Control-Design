%% TUGAS BESAR KENDALI LANJUT (KL)
% Nama  : Gede Haris Widiarta
% NIM   : 1102174038
% Title : Segway Control with LQR
% Part  : 3 (Linear Quadratic Gaussian (LQG): LQR + Kalman Filter)
%%
% Define Parameter System
km = 0.1;       %torque constant motor (Nm/A)
ke = 0.12;      %torque constant/back EMF motor (Vs/rad)
Ri = 12;        %resistor value in DC motor (ohm)
mp = 5.4;       %mass of segway robot (kg)
mw = 0.8;       %mass of segway robot wheels (kg)
l  = 1.2;       %lenght/height body form the shaft (m)
Ip = 0.048;     %moment inertia of segway robot (kg*m^2)
Iw = 0.042;     %moment inertia of segway robot wheels (kg*m^2)
g  = 9.8;       %gravity value (m/s^2)
r  = 0.15;      %radius of segway wheels (m)

B1  = (2*mw + (2*Iw/r^2) + mp);
A1  = (Ip*B1 + 2*mp*l^2*(mw + (Iw/r^2)));
%%
% Define State Space Matrix
A = [0                     1                           0            0;
     0   (2*km*ke)*(mp*l*r-Ip-mp*l^2)/Ri*r^2*A1  (mp^2*g*l^2)/A1    0;
     0                     0                           0            1;
     0      (2*km*ke)*(r*B1-mp*l)/Ri*r^2*A1      (mp*g*l*B1)/A1     0];
 
B = [0 ; 2*km*(Ip*+mp*l^2-mp*l*r)/Ri*r*A1 ; 0 ; (2*km)*(mp*l-r*B1)/Ri*r*A1];
C = [1 0 0 0];
D = [0];
%%
% LQR Control
Q = 50*eye(4);
R = 1;
K = lqr(A,B,Q,R);
%%
% Kalman Filter (Kf)
Wd = 4*eye(4);
Wn = .5;
Kf = (lqr(A',C',Wd,Wn))';
%%
% Simulink Model LGQ
t_final = 12;
sim('LQG_Sim')

t = y.time;
ytrue = y.signals.values;
yhat  = yhat.signals.values(:,1);
xtrue = x.signals.values;
xhat  = xhat.signals.values;
%%
% Plotting Response
RGB_Index = [   0      0.4470    0.7410
             0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.4940    0.1840    0.5560
             0.4660    0.6740    0.1880
             0.3010    0.7450    0.9330
             0.6350    0.0780    0.1840];

figure(1)
plot(t,ytrue,'-','LineWidth',1), hold on
plot(t,yhat,'--','LineWidth',3,'Color',RGB_Index(5,:)), grid on
title('Respon Sistem Closed-Loop dengan LQG')
xlabel('time (s)')
ylabel('x (m)')
set(legend('$y$','$\hat{y}$'),'interpreter','latex')
hold on

figure(2)
plot(t,xtrue(:,1),'-','LineWidth',1), hold on
plot(t,xhat(:,1),'--','LineWidth',3), hold on 
plot(t,xtrue(:,2),'-','LineWidth',1), hold on
plot(t,xhat(:,2),'--','LineWidth',3), hold on 
plot(t,xtrue(:,3),'-','LineWidth',1), hold on
plot(t,xhat(:,3),'--','LineWidth',3), hold on 
plot(t,xtrue(:,4),'-','LineWidth',1), hold on
plot(t,xhat(:,4),'--','LineWidth',3), hold on 
grid on
title('Respon State Closed-Loop dengan LQG')
xlabel('time (s)')
ylabel('state')
set(legend('$x$','$\hat{x}$','$\dot{x}$','$\hat{\dot{x}}$','$\theta$','$\hat{\theta}$','$\dot{\theta}$','$\hat{\dot{\theta}}$'),'Location','Northeast','interpreter','latex')