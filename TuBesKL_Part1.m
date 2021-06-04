%% TUGAS BESAR KENDALI LANJUT (KL)
% Nama  : Gede Haris Widiarta
% NIM   : 1102174038
% Title : Segway Control with LQR
% Part  : 1 (State Space, Performansi OL, Ctrb, Obsv, LQR)
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
A = [0                     1                           0              0;
     0   (2*km*ke)*(mp*l*r-Ip-mp*l^2)/Ri*r^2*A1  (mp^2*g*l^2)/A1      0;
     0                     0                           0              1;
     0      (2*km*ke)*(r*B1-mp*l)/Ri*r^2*A1      (mp*g*l*B1)/A1       0];
 
B = [0 ; 2*km*(Ip*+mp*l^2-mp*l*r)/Ri*r*A1 ; 0 ;(2*km)*(mp*l-r*B1)/Ri*r*A1];
C = [1 0 0 0;0 0 1 0];
D = [0];
%%
% Open-Loop Response
sysOL1 = ss(A,B,C(1,:),D);
sysOL2 = ss(A,B,C(2,:),D);
%PZmap
disp('Poles Sistem: ')
poles = eig(A)
subplot(2, 2, 1:4)
pzmap(sysOL1)
hold on
% pzmap(sysOL2)
grid on
title('Poles of Segway Balancing Robot System')
%Bode Plot
figure(3)
bode(sysOL1)
hold on
bode(sysOL2)
grid on
set(legend('sysOL x','sysOL \theta'))
%%
% Controllability
Matrix_Co = ctrb(A,B)
Rank_Co   = rank(ctrb(A,B))
if Rank_Co == size(Matrix_Co);
    fprintf('System is Controllable\n')
else
    fprintf('System is Uncontrollable\n')
end
%%
% Observeability
Matrix_Ob = obsv(A,C(1,:))
Rank_Ob   = rank(obsv(A,C(1,:)))
if Rank_Ob == size(Matrix_Ob);
    fprintf('System is Observeable\n')
else
    fprintf('System is Unobserveable\n')
end
%%
% Linear Quadratic Regulator (LQR)
% Define LQR Parameter
% Variasi nilai Q dan R konstan (Model 1)
for i = 1:5;
    switch i
        case 1
            Q = diag([1 1 1 1]); R = 1;
        case 2
            Q = diag([10 10 10 10]); R = 1;
        case 3
            Q = diag([30 30 30 30]); R = 1;
        case 4
            Q = diag([50 50 50 50]); R = 1;
        case 5
            Q = diag([100 100 100 100]); R = 1;
    end
    
K1 = lqr(A,B,Q,R);K1
sysCL1 = ss(A-B*K1,B,C(1,:),D);
sysCL2 = ss(A-B*K1,B,C(2,:),D);

%Plotting LQR Model 1
t = 0:0.1:20;
y1 = step(sysCL1,t);
y2 = step(sysCL2,t);

figure(4)
plot(t,y1,'-','LineWidth',1);
grid on
title('Respon Sistem Closed-Loop LQR')
xlabel('time (s)')
ylabel('position (m)')
set(legend('Q = 1','Q = 10','Q = 30','Q = 50','Q = 100'))
hold on

figure(5)
plot(t,y2,'-','LineWidth',1);
grid on
title('Respon Sistem Closed-Loop LQR')
xlabel('time (s)')
ylabel('angle (rad)')
set(legend('Q = 1','Q = 10','Q = 30','Q = 50','Q = 100'))
hold on
end

% Variasi nilai R dan Q konstan (Model 2)
for i = 1:5;
    switch i
        case 1
            Q = eye(4); R = 1;
        case 2
            Q = eye(4); R = 10;
        case 3
            Q = eye(4); R = 30;
        case 4
            Q = eye(4); R = 50;
        case 5
            Q = eye(4); R = 100;
    end
    
K2 = lqr(A,B,Q,R);K2
sysCL1 = ss(A-B*K2,B,C(1,:),D);
sysCL2 = ss(A-B*K2,B,C(2,:),D);

%Plotting LQR Model 1
t = 0:0.1:50;
y1 = step(sysCL1,t);
y2 = step(sysCL2,t);

figure(6)
plot(t,y1,'-','LineWidth',1);
grid on
title('Respon Sistem Closed-Loop LQR')
xlabel('time (s)')
ylabel('position (m)')
set(legend('R = 1','R = 10','R = 30','R = 50','R = 100'))
hold on

figure(7)
plot(t,y2,'-','LineWidth',1);
grid on
title('Respon Sistem Closed-Loop LQR')
xlabel('time (s)')
ylabel('angle (rad)')
set(legend('R = 1','R = 10','R = 30','R = 50','R = 100'))
hold on
end
        