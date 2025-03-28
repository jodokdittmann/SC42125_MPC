p.l_R = 1.59; % [m]
p.l_F = 1.1; % [m]
p.w = 1; % [m]
p.m = 1450; % [kg]
p.I_z = 2741.9; % [kg m^2]
p.R_w = 0.3; % [m]
p.h = 0.4; % [m]
p.g = 9.81; % [N kg^-1]

p.B = 7;
p.C = 1.6;
p.D = 1;

p.n_x = 6;
p.n_zeta = 3;
p.n_eta = 3;
p.n_u = 3;

p.ts = 1e-3;
p.tf = 30;

p.N = 15;

p.Q_MPC = 10*eye(p.n_zeta);
p.R_MPC = eye(p.n_u);
p.P_MPC = 100*p.Q_MPC;

p.Q_EKF = eye(p.n_x);
p.R_EKF = eye(p.n_eta);

p.n_equi = 9;

R = 10;  
n = 50;
l = 30;

p.th = 1.5;

p.X_ref = [linspace(0, l, n), R*sin(linspace(0, pi, n)) + l, l - linspace(0, l, n), - R*sin(linspace(0, pi, n))];
p.Y_ref = [zeros(1, n) + R, R*cos(linspace(0, pi, n)), zeros(1, n) - R, - R*cos(linspace(0, pi, n))];



