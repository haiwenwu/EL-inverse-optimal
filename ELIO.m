function dx = ELIO(t,x)
% % Simulation for ICCA 2019 
% % known frequency

global t_save tau_save L_save taubar_save

%% states
dx = zeros(46,1);
eta1 = x(1:6); % 6×1 
eta2 = reshape(x(7:36),6,5); eta2vec = x(7:36); % 6×5
ahat = x(37:41); % 5×1
q = x(42:43); q1 = q(1); q2 = q(2);
dq = x(44:45); dq1 = dq(1); dq2 = dq(2);
L_int = x(46);

%% System parameters
a1 = 3.9; a2 = 0.75; a3 = 1.125; b1 = 23.52; b2 = 7.35;
a = [a1; a2; a3; b1; b2];
M = [a1+2*a3*cos(q2), a2+a3*cos(q2); a2+a3*cos(q2),   a2];
C = [-a3*sin(q2)*dq2, -a3*sin(q2)*(dq1+dq2); a3*sin(q2)*dq1, 0];
g = [b1*cos(q1)+b2*cos(q1+q2); b2*cos(q1+q2)]; 

%% Design parameters
epsilon = 1;
KD = 20.0*eye(2);
lambda = 1000;
alpha = 1.0;
beta = 2;
F = [ 0     1     0     0     0     0;
      0     0     1     0     0     0;
     -1 -2.15 -1.75     0     0     0;
      0     0     0     0     1     0;
      0     0     0     0     0     1;
      0     0     0    -1 -2.15 -1.75];
H = [0 0; 0 0; 1 0; 0 0; 0 0; 0 1];
Gamma = [1 1.15 1.75    0    0    0;
         0    0    0    1 1.15 1.75]; % frequency = 1
T = [ 1.0000   -0.6101    1.3979    0    0    0;
     -0.0000   -0.3979   -0.6101    0    0    0;
      0.0000    0.6101   -0.3979    0    0    0;
      0         0         0         1.0000   -0.6101    1.3979;
      0         0         0        -0.0000   -0.3979   -0.6101;
      0         0         0         0.0000    0.6101   -0.3979];
     

%% References and disturbances
A1 = 0.3; B1 = 0.3;
sigma_1 = 1.0; sigma_2 = 1.0;
xd1 = A1 * (1 - cos(sigma_1 * t)); 
xd2 = B1 * (1 - cos(sigma_2 * t)); 
dxd1 = A1 * sigma_1 * sin(sigma_1 * t); 
dxd2 = B1 * sigma_2 * sin(sigma_2 * t); 
ddxd1 = A1 * sigma_1^2 * cos(sigma_1 * t); 
ddxd2 = B1 * sigma_2^2 * cos(sigma_2 * t); 
xd = [xd1; xd2];

e1 = q1 - xd1;
e2 = q2 - xd2;
e = [e1; e2];
de1 = dq1 - dxd1;
de2 = dq2 - dxd2;
dqr1 = dxd1 - alpha * e1;
dqr2 = dxd2 - alpha * e2;
ddqr1 = ddxd1 - alpha * de1;
ddqr2 = ddxd2 - alpha * de2;
zeta1 = dq1 - dqr1;
zeta2 = dq2 - dqr2;
zeta = [zeta1; zeta2];

sigma = 1;
d1 = sin(sigma*t) - 2 ;
d2 = sin(sigma*t) + 1 ;
d_d1 = sigma*cos(sigma*t);
d_d2 = sigma*cos(sigma*t);
dd_d1 = -sigma^2*sin(sigma*t);
dd_d2 = -sigma^2*sin(sigma*t);
d = [d1; d2];

%% functions
Yd13= (2 * ddxd1 + ddxd2)*cos(q2) - (dq2*dq1 + dq1 * dq2 + dq2 * dq2)*sin(q2);
Yd23= ddxd1 * cos(q2) + dq1 * dq1*sin(q2);
Y   = [ddxd1,       ddxd2, Yd13, cos(q1), cos(q1+q2);
          0, ddxd1+ddxd2, Yd23,       0, cos(q1+q2)];
Yce1= [0 0 sin(q2)*(-dq2*e1-(dq1-dq2)*e2) 0 0;
       0 0               sin(q2)*(dq1*e1) 0 0];
Yce2= [0, 0, sin(q2)*(-dq2*de1+dq1*de2), 0, 0;
       0, 0,    sin(q2)*((dq1-dq2)*de1), 0, 0];
Ym  = [de1,     de2, cos(q2)*(2*de1+de2), 0, 0;
        0, de1+de2,         cos(q2)*de1, 0, 0];
f1  = -(Gamma*eta2 - Y);
f0  = -(Gamma*H*Ym + alpha*Ym + alpha*Yce1);

Chat = ahat(3)*sin(q2).*[-dq2, -dq1-dq2; dq1, 0];
Mhat = [ahat(1)+2*ahat(3)*cos(q2), ahat(2)+ahat(3)*cos(q2); ahat(2)+ahat(3)*cos(q2), ahat(2)];
Psi1 = -(sqrt(alpha).*Chat - sqrt(alpha).*alpha.*Mhat - sqrt(alpha)*Gamma*H*Mhat + (1/sqrt(alpha)).*eye(2));
Psi2 = (alpha + Gamma*H)*Mhat ;

%% controller
Rni    = KD + Psi1*Psi1'/2 + Psi2'*KD^(-1)*Psi2 + Gamma*Gamma'/(2*epsilon);
tauhat = f1*ahat + Gamma*eta1;
taubar = -beta*Rni*zeta;
% Rni    = KD;
% tauhat = (f1+f0)*ahat + Gamma*eta1;
% taubar = - beta*Rni*zeta;
tau    = tauhat + taubar;

%% cost
vec_dis = [d1; d_d1; dd_d1; d2; d_d2; dd_d2];
eta_tilde = eta1 - eta2*a - T*vec_dis - H*Ym*a;
L = -2*beta*( -epsilon*(eta_tilde'*eta_tilde)- alpha* (e'*e) + e'*zeta + zeta'*Gamma*eta_tilde - zeta'*f0*ahat) + beta^2*zeta'*Rni*zeta;

%% closed-loop system
dx(1:6)   = F*eta1 + H*tau;
dx(7:36)  = reshape(F*eta2 + H*(Y+Yce2) + F*H*Ym, 30, 1);
dx(37:41) = -lambda*(f1 + f0)'*zeta;
dx(42:43) = dq;
dx(44:45) = M^(-1)*(tau - d -C*dq - g);
dx(46)    = L + taubar'*inv(Rni)*taubar;

%% save data
if(t-t_save(end) >= 0.01)
    t_save = [t_save; t];
    tau_save = [tau_save; tau'];
    L_save = [L_save; L];
    taubar_save = [taubar_save; taubar'];
end