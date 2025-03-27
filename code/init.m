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

p.n_x = 3;
p.n_zeta = 3;
p.n_u = 3;

p.ts = 1e-3;
p.tf = 30;

p.N = 15;

p.Q = 10*eye(p.n_x);
p.R = eye(p.n_u);
p.P = 10*p.Q;

p.n_equi = 15;


