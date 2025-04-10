% Car
p.l_R = 1.59; % [m]
p.l_F = 1.1; % [m]
p.w = 1; % [m]
p.m = 1450; % [kg]
p.I_z = 2741.9; % [kg m^2]
p.R_w = 0.3; % [m]
p.h = 0.4; % [m]
p.g = 9.81; % [N kg^-1]

% Tire
p.B = 7;
p.C = 1.6;
p.D = 1;

% Dimensions
p.n_x = 6;
p.n_zeta = 3;
p.n_eta = 3;
p.n_u = 3;
p.n_y = 5;
p.n_z = 3;

% Time
p.ts = 1e-2;
p.tf = 30;

% MPC
p.N = 20;
p.Q_MPC = 10*eye(p.n_zeta);
p.R_MPC = eye(p.n_u);

% KF
p.Q_KF = eye(p.n_zeta);
p.R_KF = eye(p.n_z);

% EKF
p.Q_EKF = blkdiag(0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
p.R_EKF = blkdiag(0.01, 0.01, 0.1, 0.1, 0.1);

% UKF
p.sqrtQ_UKF = chol(p.Q_EKF);
p.sqrtR_UKF = chol(p.R_EKF);
kappa = 3 - p.n_x;
alpha = 0.01;
beta = 2;
p.gamma = sqrt(alpha^2*(p.n_x + kappa));
p.Wm = 1 - p.n_x/(alpha^2*(p.n_x + kappa));
p.Wc = p.Wm + (1 - alpha^2 + beta);
p.W = 1/(2*(alpha^2*(p.n_x + kappa)));

% Path Planning
p.n_equi = 9;
p.th = 1.5;
R = 10;  
n = 50;
l = 30;
p.X_ref = [linspace(0, l, n), R*sin(linspace(0, pi, n)) + l, l - linspace(0, l, n), - R*sin(linspace(0, pi, n))];
p.Y_ref = [zeros(1, n) + R, R*cos(linspace(0, pi, n)), zeros(1, n) - R, - R*cos(linspace(0, pi, n))];



