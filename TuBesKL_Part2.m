%% TUGAS BESAR KENDALI LANJUT (KL)
% Nama  : Gede Haris Widiarta
% NIM   : 1102174038
% Title : Segway Control with LQR
% Part  : 2 (Estimator with Kalman Filter)
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
C = [1 0 0 0;0 0 1 0];
D = [0];
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
% Kalman Filter
% Define Value of Disturbances (Wd) & Noise (Wn)
Wd = eye(4); %Matrix Disturbances
Wn = 1;      %Matrix Noise

% Kalman Filter Value
Kf = (lqr(A',C(1,:)',Wd,Wn))';

% Augmented System
% xdot = Ax + Bu + Md*d + 0*n
% y    = Cx + Du + 0*d + Mn*n
Baug = [B eye(4) 0*B];                                  %[Bu Wd 0*n]
Daug = [0 0 0 0 0 1 ];                                  %[Du 0*d Wn]
sysC = ss(A,Baug,C(1,:),Daug);                          %Sistem Pengukuran
sysTruth = ss(A,Baug,eye(4),zeros(4,size(Baug,2)));     %Sistem Full Out.(dengan noise dan disturbances)
sysKF    = ss(A-Kf*C(1,:),[B Kf],eye(4),0*[B Kf]);      %Sistem Kalman Filter

% Estimasi Sistem Linear
dt = .01;
t  = dt:dt:50;

u_dist  = sqrt(Wd)*randn(4,size(t,2)); %Random Distrurbances Matrix
u_noise = sqrt(Wn)*randn(size(t));     %Random Noise Matrix
u = 0*t;
u(1/dt) = 20/dt;                       %Positive Impulse
u(15/dt) = -20/dt;                     %Negative Impulse
u_aug = [u; u_dist; u_noise];          %Input dengan Disturbance dan Noise

[y,t] = lsim(sysC,u_aug,t);            %Dengan Noise
[xtrue,t] = lsim(sysTruth,u_aug,t);    %State Sebenarnya
[xhat,t] = lsim(sysKF,[u; y'],t);      %State Estimasi

%Plotting
RGB_Index = [   0      0.4470    0.7410
             0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.4940    0.1840    0.5560
             0.4660    0.6740    0.1880
             0.3010    0.7450    0.9330
             0.6350    0.0780    0.1840];

figure(1)
plot(t,y,'Color',[.5 .5 .5])
hold on
plot(t,xtrue(:,1),'Color',[0 0 0],'LineWidth',2)
plot(t,xhat(:,1),'--','Color',RGB_Index(1,:),'LineWidth',2)
grid on
title('Estimasi Kalman Filter')
xlabel('time (s)')
ylabel('measurement')
set(legend('y (pengukuran)','y (tanpa \textit{noise})','y (estimasi dengan $\textbf{K}_f$)'),'interpreter','latex')
set(0,'defaulttextinterpreter','Latex')

figure(2)
for k=1:4
plot(t,xtrue(:,k),'-','LineWidth',1.2,'Color',RGB_Index(k,:));
hold on
plot(t,xhat(:,k),'--','LineWidth',2,'Color',RGB_Index(k,:))
end
title('Estimasi State dengan Kalman Filter (Full-Measurement)')
xlabel('time (s)')
ylabel('state')
set(legend('$x$','$\hat{x}$','$v$','$\hat{v}$','$\theta$','$\hat{\theta}$','$\omega$','$\hat{\omega}$'),'Location','NorthEast','interpreter','latex')
set(0,'defaulttextinterpreter','Latex')
grid on
