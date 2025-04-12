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
p.Q_MPC = 100*eye(p.n_zeta);
p.R_MPC = eye(p.n_u);

% KF
p.Q_KF = eye(p.n_zeta);
p.R_KF = eye(p.n_z);

% EKF
p.Q_EKF = eye(p.n_x);
p.R_EKF = eye(p.n_y);

% Path Planning
p.n_equi = 15;
p.th = 1.8;
R = 10;  
l = 60;
p.X_ref = [linspace(0, l, l), R*sin(linspace(0, pi, 3*R)) + l, l - linspace(0, l, l), - R*sin(linspace(0, pi, 3*R))];
p.Y_ref = [zeros(1, l) + R, R*cos(linspace(0, pi, 3*R)), zeros(1, l) - R, - R*cos(linspace(0, pi, 3*R))];



