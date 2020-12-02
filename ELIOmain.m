close all; clc;
clearvars -except save_t1 save_L1 save_tau1 save_taubar1 save_Jt1 save_J1

global t_save tau_save L_save taubar_save 
t_save = [0];
tau_save = [0,0];
L_save = [0];
taubar_save = [0,0];

ini = zeros(36,1);
aini = [1; 1; 1; 30; 7];
[t,x] = ode23(@ELIO,[0 20],[ini; aini; 0/180*pi; -0/180*pi; 0; 0; 0]);


lambda = 1000;
beta = 2;
A1 = 0.3; B1 = 0.3;
sigma_1 = 1.0; sigma_2 = 1.0;
xd1 = A1 * (1 - cos(sigma_1 * t)); 
xd2 = B1 * (1 - cos(sigma_2 * t)); 


figure(1)  %% error
plot(t,(x(:,42)-xd1)/pi*180,'b','LineWidth',1.5);
grid on; hold on;
plot(t,(x(:,43)-xd2)/pi*180,'--r','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Tracking error (deg)');
legend('$e_{1}$','$e_{2}$');

figure(2) % input
tau_save(1,:) = tau_save(2,:);
plot(t_save,tau_save(:,1),'b','LineWidth',1.5);
hold on
plot(t_save,tau_save(:,2),'--r','LineWidth',1.5);
grid on
xlabel('Time (s)');
ylabel('Control input (N$\cdot$m)');
legend('$\tau_{1}$','$\tau_{2}$');

figure(4) %% L
L_save(1) = L_save(2);
plot(t_save,L_save,'b','LineWidth',1.5);
grid on
xlabel('Time (s)');
ylabel('L');

figure(5) %% L_int
a1 = 3.9; a2 = 0.75; a3 = 1.125; b1 = 23.52; b2 = 7.35;
a = [a1; a2; a3; b1; b2];
cost = zeros(length(t),1);
for i=1:length(t)
   ahat = x(i,37:41);
   atilde = a' - ahat;
   cost(i) = x(i,46) + beta*(atilde*atilde')/lambda;
end
plot(t,cost,'b','LineWidth',1.5);
grid on
xlabel('Time (s)');
ylabel('Intergal of L');

figure(6) %% tau_st
taubar_save(1,:) = taubar_save(2,:);
plot(t_save,taubar_save(:,1),'b','LineWidth',1.5);
hold on; grid on;
plot(t_save,taubar_save(:,2),'--r','LineWidth',1.5);
ylabel('$\bar{\tau}$ (N$\cdot$m)');
legend('$\bar{\tau}_{1}$','$\bar{\tau}_{2}$');
xlabel('Time (s)');
